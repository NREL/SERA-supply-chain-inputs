# Set up the case to run.
using Pkg;
Pkg.add("Xpress")
Pkg.build("Xpress")
Pkg.add("Gurobi")
Pkg.build("Gurobi")

#Specify path of SERA directory

if isempty(ARGS)
    println("SERA directory not specified in input arguments, setting to default directory")
    # sera_directory = "/lustre/eaglefs/projects/sera/SERA2.0/carbon-neutral-california/SERA.jl"
    sera_directory = "/home/ccyang/SERA-supply-chain-inputs/SERA.jl"
else
    sera_directory = ARGS[1]
end


cd(sera_directory)
Pkg.activate(".");


# Set the working directory.
working_directory = @__DIR__
cd(working_directory)

# Add the SERA package to the path.
push!(LOAD_PATH, joinpath(sera_directory, "src"))

# Import modules.
using SERA,
      SERA.Supply,
      SERA.Supply.IO,
      SERA.Supply.Types,
      SERA.Supply.LookupFunctions,
      SERA.Supply.Optimize_new,
      YAML

# Read the configuation.
configuration = YAML.load(open(joinpath(working_directory, "scenario.yaml")))

# Read the network.
network = readnetwork(configuration, working_directory)
processlibrary = readprocesslibrary(configuration, working_directory)
demands = readdemands(configuration, working_directory)
periods = readperiods(configuration, working_directory)
prices = readprices(configuration, working_directory)
usage = readquantities(configuration, working_directory)
discount_rate = configuration["discountRate"]
rollingWindow = configuration["timeWindow"]
start_year = configuration["firstYear"]
last_year = configuration["lastYear"]

if in("timeWindow", keys(configuration))
    timeWindow = configuration["timeWindow"]
else
    @info "Time window not specified in scenaio.yaml file; Setting it to $(5) years"
    timeWindow = 5
end

if in("rollingStep", keys(configuration))
    rollingStep = configuration["rollingStep"]
else
    @info "Rolling step not specified in scenaio.yaml file; Setting it to $(rollingWindow) years"
    rollingStep = rollingWindow
end

if in("numberIterations", keys(configuration))
    niter = configuration["numberIterations"]
else
    @info "Number of iterations not specified in scenaio.yaml file; Setting it to $(10)"
    niter = 10
end

if in("newProduction", keys(configuration))
    build_new_production = configuration["newProduction"]
else
    @info "Building new production boolean not specified in scenaio.yaml file; Setting it to true"
    build_new_production = true
end

if in("storageAvailable", keys(configuration))
    storage_available = configuration["storageAvailable"]
else
    @info "Storage availability boolean not specified in scenaio.yaml file; Setting it to true"
    storage_available = true
end

if in("cleanProductionConstraintFiles", keys(configuration))
    clean_production_constraints = read_cleanproductionconstraints(configuration, working_directory)
else
    @info "Clean production constraints file not specified in scenaio.yaml file; SERA will not create these constraints"
    clean_production_constraints = CleanProductionConstraint[]
end

if in("productionCapacityConstraintFiles", keys(configuration))
    production_capacity_constraints = read_productioncapacityconstraints(configuration, working_directory)
else
    @info "Production capacity constraints file not specified in scenaio.yaml file; SERA will not create these constraints"
    production_capacity_constraints = ProductionCapacityConstraint[]
end

if in("storageCapacityConstraintFiles", keys(configuration))
    storage_capacity_constraints = read_storagecapacityconstraints(configuration, working_directory)
else
    @info "Storage capacity constraints file not specified in scenaio.yaml file; SERA will not create these constraints"
    storage_capacity_constraints = Dict{Year, Float64}()
end

if in("utilizationConstraintFiles", keys(configuration))
    utilization_constraints = read_utilizationconstraints(configuration, working_directory)
else
    @info "Utilization constraints file not specified in scenaio.yaml file; SERA will not create these constraints"
    utilization_constraints = Dict{Technology, Float64}()
end

if in("maximumCapacityConstraintFiles", keys(configuration))
    maximum_capacity_constraints = read_maximumcapacities(configuration, working_directory)
else
    @info "Maximum capacities file not specified in scenaio.yaml file; SERA will not create these constraints"
    maximum_capacity_constraints = Dict{Technology, Float64}()
end

if in("optimizationLog", keys(configuration))
    optimization_log = configuration["optimizationLog"]
else
    @info "Optimization logging boolean not specified in scenaio.yaml file; SERA will turn off optimization log"
    optimization_log = false
end

if in("annualization", keys(configuration))
    annualization = configuration["annualization"]
else
    @info "Annualization boolean not specified in scenaio.yaml file; SERA will turn off annualization of capital costs"
    annualization = false
end

val, t = @timed makeModel_new(network,
                            processlibrary,
                            demands,
                            prices,
                            usage,
                            periods,
                            niter,
                            discount_rate,
                            timeWindow,
                            rollingStep,
                            start_year,
                            last_year;
                            build_new_production = build_new_production,
                            storage_available = storage_available,
                            clean_production_constraints = clean_production_constraints,
                            production_capacity_constraints = production_capacity_constraints,
                            storage_capacity_constraints = storage_capacity_constraints,
                            utilization_constraints = utilization_constraints,
                            maximum_capacity_constraints = maximum_capacity_constraints,
                            optimization_log = optimization_log,
                            annualization = annualization)

printstyled("\nTHE OPTIMIZATION PROBLEM TOOK $(round(t/60, digits = 2)) MINUTES TO SOLVE. \n"; color = :green)
