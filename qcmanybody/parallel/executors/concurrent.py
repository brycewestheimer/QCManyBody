"""
Concurrent futures-based executor for parallel task execution.

This module provides an alternative executor using Python's concurrent.futures
module, which offers a simpler API than multiprocessing with both thread and
process-based execution.
"""

from typing import List, Optional, Callable
import logging
import concurrent.futures
from dataclasses import dataclass

from ..base import BaseParallelExecutor, ExecutorConfig
from ..task import ParallelTask, TaskResult

logger = logging.getLogger(__name__)


@dataclass
class ConcurrentExecutorConfig(ExecutorConfig):
    """
    Extended configuration for concurrent.futures executor.

    Attributes
    ----------
    executor_type : str
        Type of executor: "process" or "thread"
    max_workers : Optional[int]
        Maximum number of workers (None = use default from concurrent.futures)
    """

    executor_type: str = "process"
    """Type of executor: "process" (ProcessPoolExecutor) or "thread" (ThreadPoolExecutor)"""

    max_workers: Optional[int] = None
    """Maximum number of workers (None = use default from concurrent.futures)"""

    def __post_init__(self):
        """Validate configuration after initialization."""
        # Call parent validation first
        super().__post_init__()

        # Validate executor_type
        valid_types = {"process", "thread"}
        if self.executor_type not in valid_types:
            raise ValueError(
                f"Invalid executor_type '{self.executor_type}'. "
                f"Must be one of: {', '.join(sorted(valid_types))}"
            )

        # Validate max_workers
        if self.max_workers is not None and self.max_workers < 1:
            raise ValueError(f"max_workers must be >= 1 or None, got {self.max_workers}")


class ConcurrentExecutor(BaseParallelExecutor):
    """
    Executor using concurrent.futures for parallel execution.

    This executor provides an alternative to MultiprocessingExecutor using
    Python's concurrent.futures module. It supports both process-based
    and thread-based execution.

    **Process vs Thread Execution:**

    - **ProcessPoolExecutor** (default):
      - True parallelism via separate processes
      - No GIL limitations
      - Higher memory overhead
      - Best for CPU-bound tasks (QC calculations)

    - **ThreadPoolExecutor**:
      - Lightweight threads in same process
      - Subject to Python GIL
      - Lower memory overhead
      - Best for I/O-bound tasks (file reading, network)

    **Advantages over MultiprocessingExecutor:**

    - Simpler API and implementation
    - Built-in timeout handling
    - Better exception handling
    - Automatic worker cleanup

    Parameters
    ----------
    config : ConcurrentExecutorConfig
        Executor configuration

    Examples
    --------
    >>> from qcmanybody.parallel.executors import ConcurrentExecutor
    >>> from qcmanybody.parallel import ConcurrentExecutorConfig
    >>>
    >>> # Process-based execution (default)
    >>> config = ConcurrentExecutorConfig(
    ...     n_workers=4,
    ...     executor_type="process"
    ... )
    >>> executor = ConcurrentExecutor(config)
    >>>
    >>> # Thread-based execution (for I/O-bound tasks)
    >>> config = ConcurrentExecutorConfig(
    ...     n_workers=10,
    ...     executor_type="thread"
    ... )
    >>> executor = ConcurrentExecutor(config)
    """

    def __init__(self, config: Optional[ConcurrentExecutorConfig] = None):
        """
        Initialize concurrent executor.

        Parameters
        ----------
        config : ConcurrentExecutorConfig, optional
            Executor configuration. If None, uses defaults.
        """
        if config is None:
            config = ConcurrentExecutorConfig()

        # Convert ExecutorConfig to ConcurrentExecutorConfig if needed
        if not isinstance(config, ConcurrentExecutorConfig):
            config = ConcurrentExecutorConfig(
                n_workers=config.n_workers,
                timeout_per_task=config.timeout_per_task,
                max_retries=config.max_retries,
                log_level=config.log_level,
                checkpoint_interval=config.checkpoint_interval,
                checkpoint_file=config.checkpoint_file,
                cache_dir=config.cache_dir
            )

        super().__init__(config)
        self.config: ConcurrentExecutorConfig = config
        self._executor: Optional[concurrent.futures.Executor] = None
        self._n_workers: int = 0

    def initialize(self) -> None:
        """
        Initialize executor resources.

        Creates either ProcessPoolExecutor or ThreadPoolExecutor based on
        configuration.
        """
        # Determine number of workers
        if self.config.max_workers is not None:
            self._n_workers = self.config.max_workers
        elif self.config.n_workers is not None:
            self._n_workers = self.config.n_workers
        else:
            # Let concurrent.futures decide (None means use default)
            self._n_workers = 0

        # Create appropriate executor
        if self.config.executor_type == "thread":
            logger.info(f"Initializing ThreadPoolExecutor with {self._n_workers or 'default'} workers")
            self._executor = concurrent.futures.ThreadPoolExecutor(
                max_workers=self._n_workers if self._n_workers > 0 else None
            )
        elif self.config.executor_type == "process":
            logger.info(f"Initializing ProcessPoolExecutor with {self._n_workers or 'default'} workers")
            self._executor = concurrent.futures.ProcessPoolExecutor(
                max_workers=self._n_workers if self._n_workers > 0 else None
            )
        else:
            raise ValueError(
                f"Invalid executor_type '{self.config.executor_type}'. "
                "Must be 'process' or 'thread'."
            )

        self._is_initialized = True
        logger.info(f"ConcurrentExecutor initialized with {self.config.executor_type} backend")

    def shutdown(self, wait: bool = True) -> None:
        """
        Shutdown executor and clean up resources.

        Parameters
        ----------
        wait : bool
            If True, wait for all workers to finish cleanly
        """
        if self._executor is not None:
            logger.info("Shutting down concurrent.futures executor")
            self._executor.shutdown(wait=wait)
            self._executor = None
            self._is_initialized = False

    def execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None
    ) -> List[TaskResult]:
        """
        Execute tasks using concurrent.futures with dependency awareness.

        Tasks are grouped by dependency level and executed level-by-level to
        respect N-body mathematical dependencies.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to execute (pre-ordered by dependency level)
        progress_callback : Optional[Callable[[str, int, int], None]]
            Optional progress callback

        Returns
        -------
        List[TaskResult]
            Results from task execution in same order as input
        """
        if not tasks:
            return []

        if self._executor is None:
            raise RuntimeError("Executor not initialized. Use context manager or call __enter__.")

        total = len(tasks)
        logger.info(f"Executing {total} tasks with dependency awareness")

        # Import worker function
        from ..worker import execute_task

        # Group tasks by dependency level
        tasks_by_level = self._group_tasks_by_dependency_level(tasks)
        logger.info(f"Tasks organized into {len(tasks_by_level)} dependency levels")

        # Track results by task_id to preserve input order
        result_map = {}
        completed_count = 0

        timeout = self.config.timeout_per_task if self.config.timeout_per_task > 0 else None

        # Execute level by level
        for level in sorted(tasks_by_level.keys()):
            level_tasks = tasks_by_level[level]
            logger.info(f"Executing level {level} with {len(level_tasks)} tasks...")

            # Submit all tasks in this level
            future_to_task = {}
            for task in level_tasks:
                future = self._executor.submit(execute_task, task)
                future_to_task[future] = task

            # Collect results for this level as they complete
            try:
                for future in concurrent.futures.as_completed(future_to_task.keys(), timeout=timeout):
                    task = future_to_task[future]

                    try:
                        result = future.result(timeout=timeout)
                        result_map[task.task_id] = result
                        completed_count += 1

                        # Log result
                        status_symbol = "✓" if result.success else "✗"
                        logger.info(
                            f"{status_symbol} Task {completed_count}/{total} (level {level}, {task.task_id}): "
                            f"{result.execution_time:.2f}s"
                        )

                        # Progress callback
                        if progress_callback:
                            try:
                                progress_callback(task.task_id, completed_count, total)
                            except Exception as e:
                                logger.warning(f"Progress callback failed: {e}")

                        if not result.success:
                            logger.error(f"Task {task.task_id} failed: {result.error_message}")

                    except concurrent.futures.TimeoutError:
                        error_msg = f"Task {task.task_id} timed out after {timeout}s"
                        logger.error(error_msg)
                        result_map[task.task_id] = TaskResult(
                            task_id=task.task_id,
                            success=False,
                            status=TaskStatus.TIMEOUT,
                            error_message=error_msg,
                            atomic_result=None
                        )
                        completed_count += 1

                        if progress_callback:
                            try:
                                progress_callback(task.task_id, completed_count, total)
                            except Exception as e:
                                logger.warning(f"Progress callback failed: {e}")

                    except Exception as e:
                        error_msg = f"Task {task.task_id} raised exception: {str(e)}"
                        logger.error(error_msg)
                        result_map[task.task_id] = TaskResult(
                            task_id=task.task_id,
                            success=False,
                            status=TaskStatus.FAILED,
                            error_message=error_msg,
                            error_type=type(e).__name__,
                            atomic_result=None
                        )
                        completed_count += 1

                        if progress_callback:
                            try:
                                progress_callback(task.task_id, completed_count, total)
                            except Exception as e:
                                logger.warning(f"Progress callback failed: {e}")

            except concurrent.futures.TimeoutError:
                # Overall level timeout exceeded
                logger.error(f"Level {level} timeout exceeded ({timeout}s)")

                # Cancel remaining futures in this level
                for future in future_to_task.keys():
                    if not future.done():
                        future.cancel()
                        task = future_to_task[future]
                        result_map[task.task_id] = TaskResult(
                            task_id=task.task_id,
                            success=False,
                            status=TaskStatus.CANCELLED,
                            error_message="Execution cancelled due to timeout",
                            atomic_result=None
                        )
                        completed_count += 1

            logger.info(f"Level {level} complete: {len(level_tasks)} tasks finished")

        # Reconstruct results list in original task order
        results = [result_map[task.task_id] for task in tasks]

        # Summary
        successful = sum(1 for r in results if r.success)
        failed = total - successful
        logger.info(f"Execution complete: {successful} succeeded, {failed} failed out of {total} total tasks")

        return results

    def get_worker_count(self) -> int:
        """
        Get number of workers.

        Returns
        -------
        int
            Number of workers (0 if using default)
        """
        return self._n_workers


def create_concurrent_executor(
    n_workers: Optional[int] = None,
    executor_type: str = "process",
    timeout_per_task: float = 3600.0,
    max_retries: int = 2
) -> ConcurrentExecutor:
    """
    Convenience function to create a concurrent executor.

    Parameters
    ----------
    n_workers : int, optional
        Number of workers (None = use default from concurrent.futures)
    executor_type : str
        Type of executor: "process" or "thread"
    timeout_per_task : float
        Timeout per task in seconds
    max_retries : int
        Maximum number of retries for failed tasks

    Returns
    -------
    ConcurrentExecutor
        Configured executor instance

    Examples
    --------
    >>> # Process-based executor with 4 workers
    >>> executor = create_concurrent_executor(n_workers=4, executor_type="process")
    >>>
    >>> # Thread-based executor with default worker count
    >>> executor = create_concurrent_executor(executor_type="thread")
    """
    config = ConcurrentExecutorConfig(
        n_workers=n_workers,
        executor_type=executor_type,
        timeout_per_task=timeout_per_task,
        max_retries=max_retries,
        max_workers=n_workers
    )
    return ConcurrentExecutor(config)
