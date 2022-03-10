using .StandardLibrary, .BiochemicalLibrary

const DEFAULT_LIBRARY = merge(standard_library, biochemical_library)

addlibrary!(newlib) = merge!(DEFAULT_LIBRARY, newlib)