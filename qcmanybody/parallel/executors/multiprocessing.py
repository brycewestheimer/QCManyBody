"""
Multiprocessing executor for single-node parallelism.

This executor uses Python's multiprocessing module to execute tasks
in parallel across multiple CPU cores on a single machine.
"""

import multiprocessing as mp
from multiprocessing import Pool
import logging
from typing import List, Optional, Callable
import time

from ..base import BaseParallelExecutor, ExecutorConfig
from ..task import ParallelTask, TaskResult, TaskStatus
from ..worker import execute_single_task, validate_environment

logger = logging.getLogger(__name__)


def _worker_init(config: ExecutorConfig):
    """Initialize worker process.

    This function is called once when each worker process starts.
    It sets up logging and validates the environment.

    Parameters
    ----------
    config : ExecutorConfig
        Executor configuration
    """
    # Set up logging for worker
    logging.basicConfig(
        level=config.log_level,
        format="%(asctime)s - Worker %(process)d - %(levelname)s - %(message)s",
    )

    # Validate environment
    if not validate_environment():
        logger.error("Worker environment validation failed")


class MultiprocessingExecutor(BaseParallelExecutor):
    """Parallel executor using Python multiprocessing.

    This executor creates a pool of worker processes and distributes
    tasks across them. Suitable for single-node parallelism with
    CPU-bound tasks (quantum chemistry calculations).

    Parameters
    ----------
    config : Optional[ExecutorConfig]
        Executor configuration. If n_workers is None, uses all CPU cores.

    Examples
    --------
    >>> from qcmanybody.parallel import MultiprocessingExecutor, ExecutorConfig
    >>> config = ExecutorConfig(n_workers=4, timeout_per_task=1800)
    >>> executor = MultiprocessingExecutor(config)
    >>> with executor:
    ...     results = executor.execute(tasks)

    Notes
    -----
    - Each worker runs in a separate process (true parallelism)
    - Tasks must be pickleable (serializable)
    - Workers share no memory (safe from race conditions)
    - Overhead from process creation and IPC
    """

    def __init__(self, config: Optional[ExecutorConfig] = None):
        super().__init__(config)
        self._pool: Optional[Pool] = None
        self._n_workers = self.config.n_workers

    def initialize(self) -> None:
        """Create process pool and initialize workers."""
        # Auto-detect worker count if not specified
        if self._n_workers is None:
            self._n_workers = mp.cpu_count()
            self.logger.info(f"Auto-detected {self._n_workers} CPU cores")

        self.logger.info(f"Initializing MultiprocessingExecutor with {self._n_workers} workers")

        try:
            # Create process pool
            self._pool = Pool(
                processes=self._n_workers,
                initializer=_worker_init,
                initargs=(self.config,),
            )

            self._is_initialized = True
            self.logger.info("MultiprocessingExecutor initialized successfully")

        except Exception as e:
            self.logger.error(f"Failed to initialize process pool: {e}")
            raise RuntimeError(f"Executor initialization failed: {e}")

    def execute(
        self,
        tasks: List[ParallelTask],
        progress_callback: Optional[Callable[[str, int, int], None]] = None,
    ) -> List[TaskResult]:
        """Execute tasks in parallel using process pool with dependency awareness.

        Tasks are grouped by dependency level and executed level-by-level to
        respect N-body mathematical dependencies. All tasks at a given level
        execute in parallel, but the executor waits for a level to complete
        before starting the next level.

        Parameters
        ----------
        tasks : List[ParallelTask]
            Tasks to execute in parallel
        progress_callback : Optional[Callable[[str, int, int], None]]
            Progress callback function

        Returns
        -------
        List[TaskResult]
            Results in same order as input tasks

        Raises
        ------
        RuntimeError
            If executor not initialized
        """
        if not self._is_initialized or self._pool is None:
            raise RuntimeError("Executor not initialized. Call initialize() first or use context manager.")

        # Validate tasks
        self.validate_tasks(tasks)

        total = len(tasks)
        self.logger.info(f"Executing {total} tasks on {self._n_workers} workers with dependency awareness")

        start_time = time.time()

        # Group tasks by dependency level
        tasks_by_level = self._group_tasks_by_dependency_level(tasks)
        self.logger.info(f"Tasks organized into {len(tasks_by_level)} dependency levels")

        # Track results in order to preserve input task order
        result_map = {}  # task_id -> TaskResult
        completed_count = 0

        # Execute level by level
        for level in sorted(tasks_by_level.keys()):
            level_tasks = tasks_by_level[level]
            self.logger.info(f"Executing level {level} with {len(level_tasks)} tasks...")

            # Submit all tasks in this level asynchronously
            async_results = []
            for task in level_tasks:
                ar = self._pool.apply_async(
                    execute_single_task,
                    args=(task, self.config),
                    error_callback=lambda e: self.logger.error(f"Task error: {e}"),
                )
                async_results.append((task, ar))

            # Collect results for this level
            for task, ar in async_results:
                try:
                    # Wait for result with timeout
                    result = ar.get(timeout=self.config.timeout_per_task)
                    result_map[task.task_id] = result
                    completed_count += 1

                    # Log result
                    status_symbol = "✓" if result.success else "✗"
                    self.logger.info(
                        f"{status_symbol} Task {completed_count}/{total} (level {level}, {task.task_id}): "
                        f"{result.execution_time:.2f}s, worker={result.worker_id}"
                    )

                    # Progress callback
                    if progress_callback:
                        try:
                            progress_callback(task.task_id, completed_count, total)
                        except Exception as e:
                            self.logger.warning(f"Progress callback failed: {e}")

                    # Log error details if failed
                    if not result.success:
                        self.logger.error(
                            f"Task {task.task_id} failed: {result.error_type}: {result.error_message}"
                        )

                except mp.TimeoutError:
                    self.logger.error(f"Task {task.task_id} timed out after {self.config.timeout_per_task}s")
                    result_map[task.task_id] = TaskResult(
                        task_id=task.task_id,
                        success=False,
                        status=TaskStatus.TIMEOUT,
                        error_type="TimeoutError",
                        error_message=f"Task exceeded {self.config.timeout_per_task}s timeout",
                    )
                    completed_count += 1

                except Exception as e:
                    self.logger.error(f"Task {task.task_id} raised unexpected exception: {e}")
                    result_map[task.task_id] = TaskResult(
                        task_id=task.task_id,
                        success=False,
                        status=TaskStatus.FAILED,
                        error_type=type(e).__name__,
                        error_message=str(e),
                    )
                    completed_count += 1

            self.logger.info(f"Level {level} complete: {len(level_tasks)} tasks finished")

        # Reconstruct results list in original task order
        results = [result_map[task.task_id] for task in tasks]

        # Summary statistics
        total_time = time.time() - start_time
        successful = sum(1 for r in results if r.success)
        failed = total - successful
        avg_task_time = sum(r.execution_time for r in results) / total if total > 0 else 0
        speedup = sum(r.execution_time for r in results) / total_time if total_time > 0 else 1.0

        self.logger.info(
            f"Parallel execution complete: {successful} succeeded, {failed} failed\n"
            f"  Total time: {total_time:.2f}s, Avg task time: {avg_task_time:.2f}s\n"
            f"  Speedup: {speedup:.2f}x ({self._n_workers} workers)"
        )

        return results

    def shutdown(self, wait: bool = True) -> None:
        """Shutdown process pool and clean up resources.

        Parameters
        ----------
        wait : bool
            If True, wait for workers to finish cleanly.
            If False, terminate workers immediately.
        """
        if self._pool is not None:
            self.logger.info(f"Shutting down MultiprocessingExecutor (wait={wait})")

            if wait:
                # Graceful shutdown
                self._pool.close()  # No more tasks accepted
                self._pool.join()  # Wait for workers to finish
            else:
                # Forceful shutdown
                self._pool.terminate()  # Kill workers immediately
                self._pool.join()  # Clean up resources

            self._pool = None
            self._is_initialized = False
            self.logger.info("MultiprocessingExecutor shut down successfully")

    def get_info(self):
        """Get executor information."""
        info = super().get_info()
        info["n_workers"] = self._n_workers
        info["pool_active"] = self._pool is not None
        return info
