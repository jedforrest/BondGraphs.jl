# see Cobos Mendez et al. (2020), Fig. 5, Table 3, and Table 4

# Maybe 'component' should be 'element' which includes functionality for
# components and junctions. Then dispatch on the parametric type:
# - StaticStorage
# - DynamicStorage
# - Dissipator
# - Junction
# - Transformer etc.

# Doing it this way means we can get the benefit of dispatch
# while easily able to extend to new sub types for custom components
# e.g. can define a chemical energy store:
#   Ce <: StaticStorage

# These abstract types could define specific components stored in the standard library
# e.g. electric capacitor <: static storage
# these define a fixed definition that is reusable across a bond graph definition
# then the Component struct is specifically for repeated *instances* within a single bond graph model
# TODO See https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects

# TODO I should think of the model-building side of this package as generating MTK models
# from a graph description or interface - I don't need to reinvent the wheel when it comes
# to definining ports and components

abstract type BondGraphNode end

abstract type BondGraphElement <: BondGraphNode end
abstract type JunctionStructure <: BondGraphNode end

abstract type StorageElement <: BondGraphElement end
abstract type DissipatorElement <: BondGraphElement end
abstract type SourceElement <: BondGraphElement end

abstract type StaticStorageElement <: StorageElement end
abstract type DynamicStorageElement <: StorageElement end

# abstract type EffortSource <: SourceElement end
# abstract type FlowSource <: SourceElement end

abstract type ParametricJunction <: JunctionStructure end
abstract type NonParametricJunction <: JunctionStructure end

# abstract type Transformer <: ParametricJunction end
# abstract type Gyrator <: ParametricJunction end

abstract type EqualEffort <: NonParametricJunction end
abstract type EqualFlow <: NonParametricJunction end


# Julia Type trees
# using GraphRecipes, Plots
# default(size=(1000, 1000))
# plot(BondGraphVertexClass, method=:tree, fontsize=10, nodeshape=:ellipse)
