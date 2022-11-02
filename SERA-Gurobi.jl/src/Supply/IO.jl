"""
Reading and writing supply-side data.
"""
module IO


# Export functions.

export readnetwork
export readlinks
export readnodes
export readzones
export readterritories
export readexistings
export readpathways
export readprocesscosts
export readprocessinputs
export readprocesslibrary
export readprocessoutputs
export readdemands
export readdemandscenario
export readperiodfile
export readperiods
export readpricefile
export readprices
export readquantity
export readquantities
export read_cleanproductionconstraints
export read_maximumcapacities
export read_productioncapacityconstraints
export read_storagecapacityconstraints
export read_utilizationconstraints

# Import modules.

using CSV
using DataFrames
using DataFramesMeta
using SERA.Supply.Types
using SERA.Util
using Revise


"""
Read a network.

# Arguments
- `configuration`: a dictionary of configuration values
- `basedir`      : the base directory, relative to which to find data files

# Example
```
using YAML
configuration = YAML.load(open("scenario.yaml"))
network = readnetwork(configuration, "data")
```
"""
function readnetwork(configuration, basedir=".") :: Network
    files = configuration["networkFiles"]
    function readdicts(key, reader, default)
        reduce(
            merge,
            map(
                reader,
                map(
                    file -> joinpath(basedir, file),
                    files[key]
                )
            ),
            init=default
        )
    end
    Network(
        readdicts("nodeFiles", readnodes, Dict{NetworkID,Node}()),
        readdicts("linkFiles", readlinks, Dict{NetworkID,Link}()),
        readdicts("zoneFiles", readzones, Dict{Tuple{ZoneID,NetworkID},Zone}()),
        readdicts("territoryFiles", readterritories, Dict{Tuple{TerritoryID,NetworkID},Territory}()),
        readdicts("existingFiles", readexistings, Dict{Tuple{NetworkID,PeriodName},Existing}())
    )
end


"""
Read node data.

# Arguments
- `file`: the location of the node file

# Example
```
nodes = readnodes("data/nodes.tsv")
```
"""
function readnodes(file :: String) :: Dict{NetworkID,Node}
    function makeEntry(row)

        node = Node(
                                  getproperty(row, Symbol("Network ID"  )) ,
                                  getproperty(row, Symbol("X"           )) ,
                                  getproperty(row, Symbol("Y"           )) ,
                                  getproperty(row, Symbol("Area [km^2]" )) ,
            parseenum(Productive, getproperty(row, Symbol("Production?" ))),
                                  getproperty(row, Symbol("Cost [\$]"   )) ,
                                  getproperty(row, Symbol("Cost [\$/yr]"))
        )
        return (node.location => node)
    end
    Dict{NetworkID,Node}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim='\t'))
        )
    )
end


"""
Read link data.

# Arguments
- `file`: the location of the link file

# Example
```
links = readlinks("data/links.tsv")
```
"""
function readlinks(file :: String) :: Dict{NetworkID,Link}
    function makeEntry(row)
        link = Link(
                      getproperty(row, Symbol("Network ID"   )) ,
                      getproperty(row, Symbol("From Node ID" )) ,
                      getproperty(row, Symbol("To Node ID"   )) ,
                      getproperty(row, Symbol("Length [km]"  )) ,
                      getproperty(row, Symbol("Cost [\$]"    )) ,
                      getproperty(row, Symbol("Cost [\$/yr]" )) ,
                      getproperty(row, Symbol("Transmission?")),
                      getproperty(row, Symbol("Delivery?"))
        )
        return (link.location => link)
    end
    Dict{NetworkID,Link}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read zone data.

# Arguments
- `file`: the location of the zone files

# Example
```
zones = readzones("data/usa-nodes.tsv")
```
"""
function readzones(file :: String) :: Dict{Tuple{ZoneID,NetworkID},Zone}
    function makeEntry(row)
        zone = Zone(
                      getproperty(row, Symbol("Zone"       )) ,
                      getproperty(row, Symbol("Network ID" )) ,
                      getproperty(row, Symbol("Fraction"   ))
        )
        return ((zone.zone,zone.location) => zone)
    end
    Dict{Tuple{ZoneID,NetworkID},Zone}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read territory data.

# Arguments
- `file`: the location of the territory file

# Example
```
territories = readterritories("data/territories-usa.tsv")
```
"""
function readterritories(file :: String) :: Dict{Tuple{TerritoryID,NetworkID},Territory}
    function makeEntry(row)
        territory = Territory(
                      getproperty(row, Symbol("Network ID" )) ,
                      getproperty(row, Symbol("Territory"  )) ,
                      getproperty(row, Symbol("Fraction"   ))
        )
        return ((territory.territory,territory.location) => territory)
    end
    Dict{Tuple{TerritoryID,NetworkID},Territory}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read existings data.

# Arguments
- `file`: the location of the existings file

# Example
```
existings = readexistings("data/existings-SMR-scenario1.tsv")
```
"""

#function readexistings(file :: String) :: Dict{Tuple{NetworkID,PeriodName},Existing}
#    function makeEntry(row)
#        existing = Existing(
#                      getproperty(row, Symbol("Network ID"           )) ,
#                      getproperty(row, Symbol("Technology"           )) ,
#                      getproperty(row, Symbol("Year"                 )) ,
#                      getproperty(row, Symbol("Period"               )) ,
#                      getproperty(row, Symbol("Capacity [kg/yr]"     )) ,
#                      getproperty(row, Symbol("Yield [upstream/kg]"  )) ,
#                      getproperty(row, Symbol("Cost [\$/kg]"         ))
#        )
#        return ((existing.location,existing.period) => existing)
#    end
#    Dict{Tuple{NetworkID,PeriodName},Existing}(
#        map(
#            makeEntry,
#            eachrow(CSV.read(file, delim="\t"))
#        )
#    )
#end

function readexistings(file :: String) :: Dict{Tuple{NetworkID, Year, PeriodName},Existing}
    dict = Dict{Tuple{NetworkID, Year, PeriodName},Existing}()
    for row in eachrow(CSV.read(file, DataFrame; delim="\t"))
        existing = Existing(
                      getproperty(row, Symbol("Network ID"           )) ,
                      getproperty(row, Symbol("Technology"           )) ,
                      getproperty(row, Symbol("Year"                 )) ,
                      getproperty(row, Symbol("Period"               )) ,
                      getproperty(row, Symbol("Capacity [kg]"     )) ,
                      getproperty(row, Symbol("Yield [upstream/kg]"  )) ,
                      getproperty(row, Symbol("Cost [\$/kg]"         ))
        )
        existingkey = (existing.location, existing.year, existing.period)
        if haskey(dict, existingkey)
           location = existing.location
           technology = existing.technology
           year = existing.year
           period = existing.period
           capacity = dict[existingkey].capacity + existing.capacity
           yield = existing.yield
           cost = existing.cost
           dict[existingkey] = Existing(location, technology, year, period, capacity, yield, cost)
        else
           dict[existingkey] = existing
        end
    end
    return dict
end

"""
Creates Capital Cost Polynomial struct
"""
function read_capital_cost_data(row)
    column_names = names(row)

    if "Capital Cost [\$]" in column_names
        constant = getproperty(row, Symbol("Capital Cost [\$]"))
    else
        constant = 0.0
    end

    if "Capital Cost [\$/km]" in column_names
        linear = getproperty(row, Symbol("Capital Cost [\$/km]"))
    else
        linear = 0.0
    end

    if "Capital Cost [\$/km2]" in column_names
        quadratic = getproperty(row, Symbol("Capital Cost [\$/km2]"))
    else
        quadratic = 0.0
    end

    capital_cost = CapitalCostPolynomial(constant, linear, quadratic)

    return capital_cost
end


"""
Creates Storage struct
"""
function read_storage_data(row)
    column_names = names(row)
    if "Storage Technology" in column_names
        technology = getproperty(row, Symbol("Storage Technology"))
    else
        technology = "No Storage"
    end

    if "Storage [kg]" in column_names
        constant = getproperty(row, Symbol("Storage [kg]"))
    else
        constant = 0.0
    end

    if issubset(["Storage [kg/km]", "Storage [kg/km2]", "Storage [kg/km3]"], column_names)
        linear = getproperty(row, Symbol("Storage [kg/km]"))
        quadratic = getproperty(row, Symbol("Storage [kg/km2]"))
        cubic = getproperty(row, Symbol("Storage [kg/km3]"))
    else
        linear = 0.0
        quadratic = 0.0
        cubic = 0.0
    end

    storage = Storage(technology,
                      StorageCapacityPolynomial(constant, linear, quadratic, cubic))

    return storage
end

"""
Read process costs data.

# Arguments
- `file`: the location of the process costs file

# Example
```
processcosts = readprocesscosts("data/production-costs.tsv")
```
"""
function readprocesscosts(file :: String) :: Dict{ProcessKey{Float64,Nothing},ProcessCost}
    function makeEntry(row)
        processkey  = ProcessKey{Float64,Nothing}(
                         getproperty(row, Symbol("Technology"                   )),
                         getproperty(row, Symbol("Year"                         )),
                         getproperty(row, Symbol("Nameplate Capacity [kg/yr]"   )),
                         getproperty(row, Symbol("Maximum Utilization [kg/kg]"  )),
                         nothing

        )
            capital_cost = read_capital_cost_data(row)
            storage = read_storage_data(row)
            processcost = ProcessCost(
                parseenum(Productive, getproperty(row, Symbol("Production?"                         ))),
                                      getproperty(row, Symbol("Lifetime [yr]"                       )),
                                      getproperty(row, Symbol("Scaling Exponent"                    )),
                                      capital_cost,
                                      getproperty(row, Symbol("Fixed Operating Cost [fraction of CapCost/y]")),
                                      getproperty(row, Symbol("Variable Operating Cost [\$/kg]"     )),
                                      getproperty(row, Symbol("Variable Operating Cost [\$/km/kg]"  )),
                                      storage
                                      )

        return (processkey => processcost)
    end
    Dict{ProcessKey{Float64,Nothing},ProcessCost}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end

"""
Read process inputs data.

# Arguments
- `file`: the location of the process inputs file

# Example
```
processinputs = readprocessinputs("data/Inputs.tsv")
```
"""
function readprocessinputs(file :: String) :: Dict{ProcessKey{Nothing,MaterialName},ProcessInput}
    function makeEntry(row)
        processkey  = ProcessKey{Nothing,MaterialName}(
                          getproperty(row, Symbol("Technology"                   )),
                          getproperty(row, Symbol("Year"                         )),
                          getproperty(row, Symbol("Nameplate Capacity [kg/yr]"   )),
                          nothing,
                          getproperty(row, Symbol("Material"                     ))
        )
        processinput = ProcessInput(
                          getproperty(row, Symbol("Consumption [unit/kg]"     )),
                          getproperty(row, Symbol("Consumption [unit/km/kg]"  ))
        )
        return (processkey => processinput)
    end
    Dict{ProcessKey{Nothing,MaterialName},ProcessInput}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read process outputs data.

# Arguments
- `file`: the location of the process outputs file

# Example
```
processoutputs = readprocessoutputs("data/Outputs.tsv")
```
"""
function readprocessoutputs(file :: String) :: Dict{ProcessKey{Nothing,MaterialName},ProcessOutput}
    function makeEntry(row)
        processkey  = ProcessKey{Nothing,MaterialName}(
                          getproperty(row, Symbol("Technology"                   )),
                          getproperty(row, Symbol("Year"                         )),
                          getproperty(row, Symbol("Nameplate Capacity [kg/yr]"   )),
                          nothing,
                          getproperty(row, Symbol("Material"                     ))
        )
        processoutput = ProcessOutput(
                          getproperty(row, Symbol("Production [unit/kg]"     )),
                          getproperty(row, Symbol("Production [unit/km/kg]"  ))
        )
        return (processkey => processoutput)
    end
    Dict{ProcessKey{Nothing,MaterialName},ProcessOutput}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read pathway data.

# Arguments
- `file`: the location of the pathways file

# Example
```
pathways = readpathways("data/delivery-pathways.tsv")
```
"""
function readpathways(file :: String) :: Dict{PathwayKey,PathwaySpecification}
    function makeEntry(row)
        pathwaykey  = PathwayKey(
                      getproperty(row, Symbol("Pathway"              )) ,
                      getproperty(row, Symbol("Stage"                ))

        )
        pathwayspecification = PathwaySpecification(
                      getproperty(row, Symbol("Technology"           )) ,
                      getproperty(row, Symbol("Yield [upstream/kg]"  )) ,
                      getproperty(row, Symbol("Extended?"            )),
                      getproperty(row, Symbol("Transmission?"        )),
                      getproperty(row, Symbol("Delivery?"            )),
                      getproperty(row, Symbol("Format"               ))
        )
        return (pathwaykey => pathwayspecification)
    end
    Dict{PathwayKey,PathwaySpecification}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read a process library.

# Arguments
- `configuration`: a dictionary of configuration values
- `basedir`      : the base directory, relative to which to find data files

# Example
```
using YAML
configuration = YAML.load(open("scenario.yaml"))
processlibrary = readprocesslibrary(configuration, "data")
```
"""
function readprocesslibrary(configuration, basedir=".") :: ProcessLibrary
    files = configuration["processLibraryFiles"]
    function readdicts(key, reader, default)
        files_list = ["empty" for i in 1:length(files)]
        for i in 1:length(files_list)
            if key in keys(files[i])
                files_list[i] = files[i][key]
            end
        end
        filter!(x -> x !== "empty", files_list)
        reduce(
            merge,
            map(
                reader,
                map(
                    file -> joinpath(basedir, file),
                    files_list
                )
            ),
            init=default
        )
    end
    pathway_files = configuration["pathwayFiles"]
    function readarrays(reader, default)
        reduce(
            merge,
            map(
                reader,
                map(
                    file -> joinpath(basedir, file),
                    pathway_files
                )
            ),
            init=default
        )
    end
    ProcessLibrary(
        readdicts("costsFile", readprocesscosts, Dict{ProcessKey{Float64,Nothing},ProcessCost}()),
        readdicts("inputsFile", readprocessinputs, Dict{ProcessKey{Nothing,MaterialName},ProcessInput}()),
        readdicts("outputsFile", readprocessoutputs, Dict{ProcessKey{Nothing,MaterialName},ProcessOutput}()),
        readarrays(readpathways, Dict{PathwayKey,PathwaySpecification}())
    )
end

"""
Read demand data from a single file.

# Arguments
- `file`: the location of the demand file

# Example
```
demands = readdemandscenario("data/Ammonia-scenario1.tsv")
```
"""
function readdemandscenario(file :: String, demand_types::Vector)
    function makeEntry(row)
        demandkey = DemandKey(
                      getproperty(row, Symbol("Network ID"                       )) ,
                      getproperty(row, Symbol("Year"                             )) ,
                      getproperty(row, Symbol("Period"                           ))
        )

        demand = Dict(d => convert(Float64, getproperty(row, Symbol(d))) for d in demand_types)
        return (demandkey => demand)
    end
    Dict(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read in all demand files

# Arguments
- `configuration`: a dictionary of configuration values
- `basedir`      : the base directory, relative to which to find data files

# Example
```
using YAML
configuration = YAML.load(open("scenario.yaml"))
demands = readdemands(configuration, "data")
```
"""
function readdemands(configuration, basedir=".") :: Demands
    files = configuration["demandFiles"]
    dict = Dict{DemandKey,Demand}()
    for file in files
        file = joinpath(basedir, file)
        demand_types = filter!(x -> !(x in ["Network ID", "Year", "Period"]), names(CSV.read(file, DataFrame)))
        demands = readdemandscenario(file, demand_types)

        for key in keys(demands)
            if !(haskey(dict, key))
               dict[key] = demands[key]
            end
        end
    end
    return dict
end


"""
Read period data

# Arguments
- `file`: the location of the period file

# Example
```
periods = readperiodfile("data/periods.tsv")
```
"""
function readperiodfile(file :: String) :: Dict{PeriodName,Float64}
    function makeEntry(row)
        period = Period(
                      getproperty(row, Symbol("Period"            )) ,
                      getproperty(row, Symbol("Duration [yr/yr]"  ))
        )
        return (period.period => period.duration)
    end
    Dict{PeriodName,Float64}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read in all period files

# Arguments
- `configuration`: a dictionary of configuration values
- `basedir`      : the base directory, relative to which to find data files

# Example
```
using YAML
configuration = YAML.load(open("scenario.yaml"))
periods = readperiods(configuration, "data")
```
"""
function readperiods(configuration, basedir=".") :: Periods
    files = configuration["periodFiles"]
    function readarrays(reader, default)
        reduce(
            merge,
            map(
                reader,
                map(
                    file -> joinpath(basedir, file),
                    files
                )
            ),
            init=default
        )
    end
    Periods(
        readarrays(readperiodfile, Dict{PeriodName,Float64}())
    )
end


"""
Read price data from a single file.

# Arguments
- `file`: the location of the price file

# Example
```
prices = readpricefile("data/salt-cavern-non-billable.tsv")
```
"""
function readpricefile(file :: String) :: Dict{PriceKey,PriceCurve}
    dict = Dict{PriceKey,PriceCurve}()
    for row in eachrow(CSV.read(file, DataFrame; delim="\t"))
        pricekey = PriceKey(
                      getproperty(row, Symbol("Material" )) ,
                      getproperty(row, Symbol("Zone"     )) ,
                      getproperty(row, Symbol("Year"     )) ,
                      getproperty(row, Symbol("Period"     ))
        )
        pricedata = PriceData(
                      getproperty(row, Symbol("Price [\$/unit]"      )) ,
                      getproperty(row, Symbol("Billable?"  ))
        )
        quantity = getproperty(row, Symbol("Quantity [unit]"))
        if haskey(dict, pricekey)
           pricecurve = dict[pricekey]
           pricecurve[quantity] = pricedata
        else
           pricecurve = PriceCurve(getproperty(row, Symbol("Quantity [unit]")) => pricedata)
        end
        dict[pricekey] = pricecurve
    end
    return dict
end


"""
Read material usage from a single price file.

# Arguments
- `file`: the location of the price file

# Example
```
materialusage = readquantity("data/salt-cavern-non-billable.tsv")
```
"""
function readquantity(file :: String) :: Dict{PriceKey,Float64}
    function makeEntry(row)
        pricekey = PriceKey(
                      getproperty(row, Symbol("Material"  )) ,
                      getproperty(row, Symbol("Zone"      )) ,
                      getproperty(row, Symbol("Year"      )) ,
                      getproperty(row, Symbol("Period"      ))
        )
        return (pricekey => getproperty(row, Symbol("Quantity [unit]")))
    end
    Dict{PriceKey,Float64}(
        map(
            makeEntry,
            eachrow(CSV.read(file, DataFrame; delim="\t"))
        )
    )
end


"""
Read in all price files

# Arguments
- `configuration`: a dictionary of configuration values
- `basedir`      : the base directory, relative to which to find data files

# Example
```
using YAML
configuration = YAML.load(open("scenario.yaml"))
prices = readprices(configuration, "data")
```
"""
function readprices(configuration, basedir=".") :: Prices
    files = configuration["priceFiles"]
    function readarrays(reader, default)
        reduce(
            merge,
            map(
                reader,
                map(
                    file -> joinpath(basedir, file),
                    files
                )
            ),
            init=default
        )
    end
    Prices(
        readarrays(readpricefile, Dict{PriceKey,PriceCurve}())
    )
end

"""
Read in quantities from all price files

# Arguments
- `configuration`: a dictionary of configuration values
- `basedir`      : the base directory, relative to which to find data files

# Example
```
using YAML
configuration = YAML.load(open("scenario.yaml"))
usage = readquantities(configuration, "data")
```
"""
function readquantities(configuration, basedir=".") :: MaterialUsage
    files = configuration["priceFiles"]
    function readarrays(reader, default)
        reduce(
            merge,
            map(
                reader,
                map(
                    file -> joinpath(basedir, file),
                    files
                )
            ),
            init=default
        )
    end
    MaterialUsage(
        readarrays(readquantity, Dict{PriceKey,Float64}())
    )
end

function read_cleanproductionconstraints(configuration, basedir=".")
    clean_production_constraints = CleanProductionConstraint[]
    if in("cleanProductionConstraintFiles", keys(configuration))
        files = configuration["cleanProductionConstraintFiles"]
        for file in files
            file_path = joinpath(basedir, file)
            for row in eachrow(CSV.read(file_path, DataFrame; delim="\t"))
                push!(clean_production_constraints,
                    CleanProductionConstraint(row["Year"] ,
                                        split(row["Eligible Technologies"], ", ") ,
                                        row["Percentage Requirement"])
                    )
            end
        end

    end

    return clean_production_constraints
end

function read_productioncapacityconstraints(configuration, basedir=".")
    capacity_constraints = ProductionCapacityConstraint[]
    if in("productionCapacityConstraintFiles", keys(configuration))
        files = configuration["productionCapacityConstraintFiles"]

        for file in files
            file_path = joinpath(basedir, file)
            for row in eachrow(CSV.read(file_path, DataFrame; delim="\t"))
                push!(capacity_constraints,
                        ProductionCapacityConstraint( row["Network ID"],
                                        row["Year"] ,
                                        split(row["Eligible Technologies"], ", ") ,
                                        row["Capacity [kg/year]"]
                                        )
                    )
            end
        end

    end

    return capacity_constraints
end

function read_storagecapacityconstraints(configuration, basedir=".")
    capacity_constraints = Dict{Year, Float64}()
    if in("storageCapacityConstraintFiles", keys(configuration))
        files = configuration["storageCapacityConstraintFiles"]

        for file in files
            file_path = joinpath(basedir, file)
            for row in eachrow(CSV.read(file_path, DataFrame; delim="\t"))
                capacity_constraints[row["Year"]] = row["Capacity [kg]"]
            end
        end
    end

    return capacity_constraints
end

function read_utilizationconstraints(configuration, basedir=".")
    utilization_constraints = Dict{Technology, Float64}()
    if in("utilizationConstraintFiles", keys(configuration))
        files = configuration["utilizationConstraintFiles"]

        for file in files
            file_path = joinpath(basedir, file)
            for row in eachrow(CSV.read(file_path, DataFrame; delim="\t"))
                utilization_constraints[row["Technology"]] = row["Utilization Requirement"]
            end
        end

    end

    return utilization_constraints
end

function read_maximumcapacities(configuration, basedir=".")
    maximum_capacity_constraints = Dict{Technology, Float64}()
    if in("maximumCapacityConstraintFiles", keys(configuration))
        files = configuration["maximumCapacityConstraintFiles"]

        for file in files
            file_path = joinpath(basedir, file)
            for row in eachrow(CSV.read(file_path, DataFrame; delim="\t"))
                maximum_capacity_constraints[row["Technology"]] = row["Maximum Capacity [kg/yr]"]
            end
        end

    end

    return maximum_capacity_constraints
end
end
