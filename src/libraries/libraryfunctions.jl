using .StandardLibrary, .BiochemicalLibrary

const DEFAULT_LIBRARY = merge(standard_library, biochemical_library)

"""
    addlibrary!(newlib)

Combine the library `newlib` with the default library used within BondGraphs. `newlib` will
need to be in the form of a dictionary, and new components should follow the below schema.

NOTE: This library is likely to change in the future. Do not rely on this schema too much.

## Library Schema
- description -> written description of the component and definitions
- numports -> the number of ports in the component
- variables ->
  - parameters -> constant parameters, unique for each component instance
  - globals -> global parameters (i.e. no namespace)
  - states -> time-dependent state variables
  - controls -> time-dependent parameters that can accept julia functions
- equations -> symbolic description of the constitutive equations
"""
addlibrary!(newlib) = merge!(DEFAULT_LIBRARY, newlib)

# TODO: docstrings in the libraries should remove the need for this function
"""
    description(comp::Symbol)
    description(lib, comp::Symbol)

Print the description of the component with symbol `comp` in library `lib`.
"""
description(comp::Symbol) = description(DEFAULT_LIBRARY, comp)
description(lib, comp::Symbol) = print(lib[comp][:description])
