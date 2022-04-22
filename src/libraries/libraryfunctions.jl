using .StandardLibrary, .BiochemicalLibrary

const DEFAULT_LIBRARY = merge(standard_library, biochemical_library)

addlibrary!(newlib) = merge!(DEFAULT_LIBRARY, newlib)

description(comp::Symbol) = description(DEFAULT_LIBRARY, comp)
description(lib, comp::Symbol) = print(lib[comp][:description])