using .StandardLibrary

const DEFAULT_LIBRARY = Ref(standard_library)

set_library!(lib = standard_library) = DEFAULT_LIBRARY[] = lib