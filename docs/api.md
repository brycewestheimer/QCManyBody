::: qcmanybody.ManyBodyCore
    options:
        show_root_heading: true

::: qcmanybody.ManyBodyComputer
    options:
        show_root_heading: true
        inherited_members: true
        members:
          - from_manybodyinput

$pydantic: qcmanybody.computer.ManyBodyComputer

::: qcmanybody.utils
    options:
        show_root_heading: true

::: qcmanybody.builder
    options:
        show_root_heading: true

::: qcmanybody.models.hierarchy
    options:
        show_root_heading: true
        members:
          - FragmentHierarchy
          - HMBESpecification
          - SchengenSpecification

::: qcmanybody.hmbe_filter
    options:
        show_root_heading: true
        members:
          - passes_hmbe_filter
          - filter_compute_list
          - select_schengen_terms
          - compute_distance_metric

::: qcmanybody.hmbe_enumerate
    options:
        show_root_heading: true
        members:
          - enumerate_hmbe_terms
