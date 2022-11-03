"""
Supply-side datatypes.
"""
module Types


# Export datatypes.

using Base: Float64
export CapitalCostPolynomial
export CleanProductionConstraint
export Construction
export Demand
export DemandKey
export Demands
export Existing
export InfrastructureID
export Link
export MaterialName
export MaterialUsage
export Network
export NetworkID
export Node
export Pathway
export PathwayKey
export PathwaySpecification
export Period
export Periods
export PeriodName
export PriceCurve
export PriceData
export PriceKey
export Prices
export ProcessCost
export ProcessInput
export ProcessKey
export ProcessLibrary
export ProcessOutput
export ProductionCapacityConstraint
export Productive
export Storage
export StorageCapacityPolynomial
export Technology
export Territory
export TerritoryID
export Year
export Zone
export ZoneID


"""
Identifier for network locations.
"""
const NetworkID = String

"""
Identifier for network zones
"""
const ZoneID = String

"""
Identifier for network territories
"""
const TerritoryID = String

"""
Identifier for year
"""
const Year = Int64

"""
Name for Production Technology/Infrastructure Component
"""
const Technology = String


"""
Name for a material.
"""
const MaterialName = String

"""
Name for Subannual Period
"""
const PeriodName = String

"""
Identifier for infrastructure
"""
const InfrastructureID = String

"""
Constraints on production.
"""
@enum Productive Yes Onsite Central No None


"""
Network nodes.
"""
struct Node
    "identifier for the node, which must be unique among nodes and links"
    location   :: NetworkID
    "x-coordinate for the node"
    x          :: Float64
    "y-coordinate for the node"
    y          :: Float64
    "area of the demand center, or zero for other locations"
    area       :: Float64
    "whether and what type of production technology can be placed at the node"
    productive :: Productive
    "cost of purchasing the land"
    sale       :: Float64
    "cost of renting the land"
    rent       :: Float64
end


"""
Network links.
"""
struct Link
    "identifier for the link, which must be unique among nodes and links"
    location     :: NetworkID
    "unique identifier for the node at the start of the link"
    from         :: NetworkID
    "unique identifier for the node at the end of the link"
    to           :: NetworkID
    "the length of the link"
    length       :: Float64
    "cost of purchasing the land"
    sale         :: Float64
    "cost of renting the land"
    rent         :: Float64
    "whether this link can be used for long distance transmission of hydrogen"
    transmission :: Bool
    "whether this link can be used for short distance delivery of hydrogen"
    delivery     :: Bool
end


"""
Network zones.
"""
struct Zone
    "geographic zone for the price or intensity"
    zone         :: ZoneID
    "identifier for the node or link"
    location     :: NetworkID
    "Fraction of the node or link residing in the zone"
    fraction     :: Float64
end

"""
Network territories.
"""
struct Territory
    "identifier for infrastructure in the service territory"
    location   :: NetworkID
    "unique identifier for the service territory"
    territory  :: TerritoryID
    "Fraction of the node or link residing in the territory"
    fraction   :: Float64
end


"""
Network existings.
"""
struct Existing
    "the location of the existing infrastructure"
    location    :: NetworkID
    "name of the production technology"
    technology  :: Technology
    "the year when the infrastructure becomes available"
    year        :: Year
    "subannual period, ordered lexicographically"
    period      :: PeriodName
    "production, transmission, delivery, or dispensing capacity"
    capacity    :: Float64
    "the output of hydrogen per unit input of hydrogen (e.g., losses), only relevant for links"
    yield       :: Float64
    "the cost of producing hydrogen at the node or transporting it across the link"
    cost        :: Float64
end


"""
Networks
"""
struct Network
    "Dictionary of nodes."
    nodes :: Dict{NetworkID,Node}
    "Dictionary of links."
    links :: Dict{NetworkID,Link}
    "Dictionary of zones."
    zones :: Dict{Tuple{ZoneID,NetworkID},Zone}
    "Dictionary of territories."
    territories :: Dict{Tuple{TerritoryID,NetworkID},Territory}
    "Dictionary of existings."
    existings :: Dict{Tuple{NetworkID, Year, PeriodName},Existing}
end


"""
Name for Pathway
"""
const Pathway = String

"""
Process Key
"""
struct ProcessKey{D,M}
    "name of the production technology"
    technology :: Technology
    "time period for the cost data"
    year       :: Year
    "production, transmission, delivery, or dispensing capacity"
    nameplate  :: Float64
    "the maximum fraction of the capacity that can be utilized over the time period"
    dutycycle  :: D
    "the material"
    material   :: M
end

"""
Capital Cost Polynomial Coefficients
"""
struct CapitalCostPolynomial
    "constant term in the polynomial"
    constant :: Float64
    "linear term (w.r.t to length) in the polynomial"
    linear :: Float64
    "quadratic term (w.r.t to length) in the polynomial"
    quadratic :: Float64
end

"""
Storage Capacity Polynomial Coefficients
"""
struct StorageCapacityPolynomial
    "constant term in the polynomial"
    constant :: Float64
    "linear term (w.r.t to length) in the polynomial"
    linear :: Float64
    "quadratic term (w.r.t to length) in the polynomial"
    quadratic :: Float64
    "cubic term (w.r.t to length) in the polynomial"
    cubic :: Float64
end

"""
Storage
"""
struct Storage
    "storage technology"
    technology  ::  Technology
    "storage capacity polynomial"
    capacitypolynomial  :: StorageCapacityPolynomial
end


"""
Process Cost
"""
struct ProcessCost
    "whether and what type of production technology is present"
    productive          :: Productive
    "the working lifetime of the technology"
    lifetime            :: Int64
    "the exponent for economies of scale"
    scaling             :: Float64
    "capital cost polynomial coefficients"
    capitalcost         :: CapitalCostPolynomial
    "fixed cost as a fraction of capital cost"
    fixedcostfrac           :: Float64
    "intercept of the variable operating cost equation"
    variablecost        :: Float64
    "coefficient of distance in the variable operating cost equation"
    variablecoststretch :: Float64
    "the storage technology, if any, which must be also be separately listed here"
    storage             :: Union{Storage, Nothing}
end


"""
Process Input
"""
struct ProcessInput
    "amount of material consumed, per kg of hydrogen output of the process"
    consumptionrate        :: Float64
    "coefficient of distance in amount of material consumed"
    consumptionratestretch :: Float64
end


"""
Process Output
"""
struct ProcessOutput
    "amount of material produced, per kg of hydrogen output of the process"
    productionrate        :: Float64
    "coefficient of distance in amount of material produced"
    productionratestretch :: Float64
end


"""
Pathway Key
"""
struct PathwayKey
    "name of the pathway"
    pathway :: Pathway
    "sequence, counting from upstream to downstream of the component in the pathway"
    stage   :: Int8
end


"""
Pathway Specification
"""
struct PathwaySpecification
    "name of the infrastructure component"
    technology   :: Technology
    "the output of hydrogen per unit input of hydrogen (due to losses)"
    yield        :: Float64
    "whether the component is a link extending between nodes"
    extended     :: Bool
    "whether this stage is for long distance transmission of hydrogen"
    transmission :: Bool
    "whether this stage is for short distance delivery of hydrogen"
    delivery     :: Bool
    "the form of hydrogen, typically GH2 or LH2"
    format       :: String
end


"""
Process descriptions.
"""
struct ProcessLibrary
    "Dictionary of process costs"
    costs    :: Dict{ProcessKey{Float64,Nothing},ProcessCost}
    "Dictionary of process inputs"
    inputs   :: Dict{ProcessKey{Nothing,MaterialName},ProcessInput}
    "Dictionary of process outputs"
    outputs  :: Dict{ProcessKey{Nothing,MaterialName},ProcessOutput}
    "Dictionary of pathways"
    pathways :: Dict{PathwayKey,PathwaySpecification}
end


"""
Demand Key.
"""
struct DemandKey
    "the location of the demand"
    location           :: NetworkID
    "the year in question"
    year               :: Year
    "subannual period, ordered lexicographically"
    period             :: PeriodName
end


"""
Demand.
"""
const Demand = Dict{String, Float64}

const Demands = Dict{DemandKey,Demand}


"""
Periods.
"""
struct Period
    "subannual period, ordered lexicographically"
    period   :: PeriodName
    "fraction of the year residing in this period"
    duration :: Float64
end


const Periods = Dict{PeriodName,Float64}


"""
Key for prices.
"""
struct PriceKey
    "name of the feedstock, including its unit of measure"
    material :: MaterialName
    "geographic zone for the price"
    zone     :: ZoneID
    "year for the price"
    year     :: Year
    "intra-year period for the price"
    period  :: PeriodName
end


"""
Prices in price curves.
"""
struct PriceData
    "price in dollars per unit"
    price    :: Float64
    "whether the price is to be included in hydrogen pricing"
    billable :: Bool
end


"""
Supply curve for prices.
"""
const PriceCurve = Dict{Float64,PriceData}


const Prices = Dict{PriceKey,PriceCurve}


"""
Record of how much material is used in each zone in each year.
"""
const MaterialUsage = Dict{PriceKey,Float64}


"""
Construction
"""
mutable struct Construction
    "unique name of the infrastructure"
    infrastructureid :: InfrastructureID
    "the location of the infrastructure"
    networkid        :: NetworkID
    "name of the production technology"
    technology :: Technology
    "whether and what type of production technology can be placed at the node"
    productive          :: Productive
    "the year the infrastructure was constructed"
    year       :: Year
    "the working lifetime of the technology"
    lifetime            :: Int64
    "production, transmission, delivery, or dispensing capacity"
    nameplate  :: Float64
    "the maximum fraction of the capacity that can be utilized over the time period"
    dutycycle  :: Float64
    "the length of the infrastructure"
    length :: Float64
    "the capital cost of the infrastructure"
    capitalcost         :: Float64
    "the fixed operating cost of the infrastructure"
    fixedcost           :: Float64
    "the variable operating cost of the infrastructure"
    variablecost        :: Float64
    "the total available storage capacity"
    storagecapacity     :: Float64
end


"""
Flow
"""
struct Flow
    "unique name of the infrastructure"
    infrastructureid :: InfrastructureID
    "the year in question"
    year             :: Year
    "the amount of hydrogen produced by the infrastructure"
    production       :: Float64
    "the amount of hydrogen flowing through the infrastucture"
    flow             :: Float64
    "the amount of hydrogen lost by the infrastructure"
    loss             :: Float64
    "the total cost associated with the infrastructure for the year in question"
    cost             :: Float64
    "the remaining capital value if the infrastructure where salvaged"
    salvage          :: Float64
end

struct CleanProductionConstraint
    "the year in question"
    year            :: Year
    "technologies eligible for participation"
    eligible_technologies :: Vector{Technology}
    "minimum percentage of H2 demand to be met by eligible technologies"
    requirement           :: Float64
end

struct ProductionCapacityConstraint
    "the location of production capacity"
    location           :: NetworkID
    "the year in question"
    year            :: Year
    "technologies eligible for participation"
    eligible_technologies :: Vector{Technology}
    "minimum accumulated capacity requirement to be met by eligible technologies"
    capacity           :: Float64
end

end
