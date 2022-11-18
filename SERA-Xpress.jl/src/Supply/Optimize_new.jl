module Optimize_new

using Base: String, Float64
using CSV
using Distributed
using DataFrames
using DataFramesMeta
using SERA.Supply.Types
using SERA.Supply.LookupFunctions
using SERA.Util
using JuMP
using GLPK
using Ipopt
using Xpress
using Gurobi
using AxisArrays

import MathOptInterface
import MathOptInterface: AbstractOptimizer

const MOI = MathOptInterface

import SERA.Supply.Types: MaterialName, MaterialUsage, Prices, ZoneID, Year, Network, ProcessLibrary, Demand, No, None, PeriodName, PathwayKey, Pathway, Construction

export makeModel_new

include("ObjectiveFunction.jl")

# Function that yields next period (season) in a given year or the first period of the next year

function next((year,period), period_set, year_set)
  if ((length(period_set) == 1) || (period == period_set[end]))
     next_year = year_set[findfirst(x -> x == year, year_set) + 1]
     result  = (next_year, period_set[1])
  else
     ind = findfirst(isequal(period), period_set)
     result = (year, period_set[ind+1])
  end
  return(result)
end

# Function that yields previous period (season) in a given year or the final period of the previous year

function prev((year,period), period_set, year_set)
  if ((length(period_set) == 1) || (period == period_set[1]))
     prev_year = year_set[findfirst(x -> x == year, year_set) - 1]
     result  = (prev_year, period_set[end])
  else
     ind = findfirst(isequal(period), period_set)
     result = (year, period_set[ind-1])
  end
  return(result)
end

function current_and_prev_years(year, year_set)
    result = filter(y -> y <= year, year_set)
    return(result)
  end

function prev_years(year, year_set)
    result = filter(y -> y < year, year_set)

    if isempty(result)
        error("no year previous to $(year)")
    end
    return(result)
end

# Function to extract subset of data frame by variable, node, year and period
function subdf0(flowDF, variable, node, year, period)
    sub0 = (flowDF.Variable .== variable) .&  (flowDF.Node .== node) .& (flowDF.Year .== year) .& (flowDF.Period .== period)
    return(sub0)
  end

# Function to extract subset of data frame by variable, technology, node, and year
function subdf1(flowDF, variable, technology, node, year)
    sub1 = (flowDF.Variable .== variable) .&  (flowDF.Technology .== technology) .& (flowDF.Node .== node) .& (flowDF.Year .== year)
    return(sub1)
end

# Function to extract subset of data frame by pathway and stage
function subdf2(flowDF, pathway, stage)
    sub2 = (flowDF.Pathway .== pathway) .&  (flowDF.Stage .== stage)
    return(sub2)
end

# Function to extract subset of data frame by period and link
function subdf3(flowDF, period, link)
    sub3 = (flowDF.Period .== period) .&  (flowDF.Link .== link)
    return(sub3)
end

# Function to extract subset of data frame by period and year
function subdf4(flowDF, year_period)
    sub4 = (flowDF.Period .== year_period[2]) .&  (flowDF.Year .== year_period[1])
    return(sub4)
end

# Function to extract subset of data frame byvariable
function subdf5(flowDF, variable)
    sub5 = (flowDF.Variable .== variable)
    return(sub5)
end


# Function to extract subset of data frame by infrastructure id
function subdf_id(flowDF, infrastructure_id)
    sub = (flowDF.Infrastructure_ID .== infrastructure_id)
    return(sub)
end

# Formulate Optimization Model

function makeModel_new(network::Network,
                       processlibrary::ProcessLibrary,
                       demands::Demands, prices::Prices,
                       usage::MaterialUsage, periods::Periods,
                       niter::Int64, discount_rate::Float64,
                       rolling_window::Int64,
                       rolling_step::Int64,
                       start_year::Int64,
                       last_year::Int64;
                       build_new_production = true,
                       storage_available::Bool = true,
                       clean_production_constraints::Vector{CleanProductionConstraint} = CleanProductionConstraint[],
                       production_capacity_constraints::Vector{ProductionCapacityConstraint} = ProductionCapacityConstraint[],
                       production_capacity_switch::String = "Nodal",
                       storage_capacity_constraints::Dict{Year, Float64} = Dict{Year, Float64}(),
                       utilization_constraints::Dict{Technology, Float64} = Dict{Technology, Float64}(),
                       maximum_capacity_constraints::Dict{Technology, Float64} = Dict{Technology, Float64}(),
                       optimization_log::Bool = true,
                       annualization::Bool = false,
                       )

    ϵ_max = 0.001
    min_pipeline_distance = 20

    @assert rolling_step <= rolling_window

    Xpress_optimizer = optimizer_with_attributes(Xpress.Optimizer,
                                                "MIPRELSTOP" => 1e-10,
                                                "BARGAPSTOP" => 1e-10,
                                                "BARDUALSTOP" => 1e-10,
                                                "BARPRIMALSTOP" => 1e-10,
                                                "MATRIXTOL" => 1e-10,
                                                "BARORDER" => 1,
                                                "OUTPUTLOG" => optimization_log,
                                                "DEFAULTALG" => 4,
                                                "CROSSOVER"=>0)
    Gurobi_optimizer = optimizer_with_attributes(Gurobi.Optimizer,
                                                "OutputFlag" => optimization_log,
                                                "BarConvTol" => 1e-9,
                                                "BarOrder" => 1,
                                                "BarHomogeneous" => 1, # set to one to detect infeasibilities and unboundedness
                                                "NumericFocus" => 3, # Set to the maximum setting to try avoiding numerical issues
                                                "Method" => 2,   # barrier method
                                                "Crossover" => 0)       # -1 means it chooses cross over automoatically

    printstyled("STARTING DATA PRE-PROCESSING \n"; color = :yellow)
    #model = Model(with_optimizer(Ipopt.Optimizer, max_iter=100))
    ################################## YEAR AND PERIOD SETS ##################################################

    println("Processing Temporal Data")
    # Create year, period and (year, period) sets
    demands_year_set = Vector{Year}()
    demands_period_set = Vector{PeriodName}()
    demands_year_period_set = Vector{Tuple{Year,PeriodName}}()

    for key in keys(demands)
	    push!(demands_year_set, key.year)
	    push!(demands_period_set, key.period)
        push!(demands_year_period_set, (key.year, key.period))
    end

    total_year_set = sort(unique(demands_year_set))
    period_set = sort(unique(demands_period_set))
    total_year_period_set = sort(unique(vcat(demands_year_period_set)))

    filter!(x -> (x >= start_year) && (x <= last_year), total_year_set)
    filter!(x -> (x[1] >= start_year) && (x[1] <= last_year), total_year_period_set)


    total_year_set_string = string.(total_year_set)

    npv_array = Dict(y => 1/(1 + discount_rate) ^ (y - start_year) for y in total_year_set)


    min_period_interval = minimum(sort(unique(values(periods))))

    ################################## NODES ##################################################

    println("Processing Nodes")
    # Node numbers
    # All nodes in the network
    node_data = network.nodes
    node_set = collect(keys(node_data))

    # All production technologies in the network
    selected_costs = filter(kv->!(kv[2].productive in [No, None]), processlibrary.costs)
    production_set = sort(unique(getfield.(keys(selected_costs), :technology)))

    # All existing technologies in the network
    existings = network.existings
    existings_set = sort(unique(getfield.(values(existings), :technology)))

    # All nodes with existing capacity
    existings_node_set = sort(unique(getfield.(values(existings), :location)))

    # All nodes with demands
    demands_node_set = sort(unique(getfield.(keys(demands), :location)))

    #=
    # Reset node set to only those nodes with demands and existings
    if (length(unique(vcat(demands_node_set, existings_node_set))) < length(node_set))
        node_set = unique(vcat(demands_node_set, existings_node_set))
    end
    =#

    # Dictionary connecting nodes to zones
    dict_zones = Dict{}()
    node_existings = Dict(n => String[] for n in node_set)

    # Create a dictionary that maps each node to corresponding zones and the fraction of the node in each zone
    for node in node_set
        dict_node_zones = filter(kv ->(kv[2].location == node), network.zones)
        node_zones = [(x.zone, x.fraction) for x in values(dict_node_zones)]
        dict_zones[node] = node_zones

        for existing in values(existings)
            if getfield(existing, :location) == node
                tech_name = getfield(existing, :technology)
                if !(tech_name in node_existings[node])
                    push!(node_existings[node], tech_name)
                end
            end
        end
    end

    ################################## DEMANDS ###############################################

    println("Processing Demand")
    demands_total = Dict{Tuple{NetworkID,Tuple{Year,PeriodName}},Float64}()

    for (key, value) in demands
        demands_total[(key.location, (key.year, key.period))] = sum(values(value))
    end

    problematic_nodes = String[]
    #=
    for n in demands_node_set
        for y in total_year_set
            for t in period_set
                if !in((n, (y, t)), keys(demands_total))
                    demands_total[(n, (y, t))] = 0.0
                end

                if y > start_year
                    years = collect(start_year:y - 1)
                    for x in years
                        if demands_total[(n, (y, t))] < 0.1 * demands_total[(n, (x, t))]
                            if !(n in problematic_nodes)
                                    push!(problematic_nodes, n)
                            end
                        end
                    end
                end
            end
        end
    end
    println(problematic_nodes)
    quit()
    =#

    for n in demands_node_set
        for y in total_year_set
            for t in period_set
                if !in((n, (y, t)), keys(demands_total))
                    demands_total[(n, (y, t))] = 0.0
                end
            end
        end
    end

    ############################ PRODUCTION CAPABILITY ########################################
    production_nodes = filter(n -> node_data[n].productive != No && node_data[n].productive != None, node_set)

    productive_options = Dict{Productive, Vector{Technology}}()
    node_production_set = Dict{String, Vector{Technology}}()
    onsite_nodes = Vector{String}()

    for m in instances(Productive)
        productive_options[m] = filter(t -> occursin(string(m), t), production_set)
    end

    for n in production_nodes
        productive = node_data[n].productive
        node_production_set[n] = productive_options[productive]
        if string(productive) == "Onsite"
            push!(onsite_nodes, n)
            @assert n in demands_node_set
        end
    end

    ################################## EXISTINGS #############################################
    println("Processing Existings")
    existings_capacity = Dict{Tuple{NetworkID, Technology, Tuple{Year,PeriodName}}, Float64}()
    existings_yield = Dict{Tuple{NetworkID, Technology, Tuple{Year,PeriodName}}, Float64}()
    existings_cost = Dict{Tuple{NetworkID, Technology, Tuple{Year,PeriodName}}, Float64}()

    for value in values(existings)
        existings_capacity[(value.location, value.technology, (value.year, value.period))] = value.capacity
        existings_yield[(value.location, value.technology, (value.year, value.period))] = value.yield
        existings_cost[(value.location, value.technology, (value.year, value.period))] = value.cost
    end

    existings_keys  = keys(existings_capacity)

    for n in node_set
        for x in node_existings[n]
            for t in total_year_set
                for p in period_set
                    if !((n, x, (t, p)) in existings_keys)
                        existings_capacity[(n, x, (t, p))] = 0.0
                        existings_yield[(n, x, (t, p))] = 0.0
                        existings_cost[(n, x, (t, p))] = 0.0
                    end
                end
            end
        end
    end

    ################################# LINKS ##################################################
    println("Processing Links")
    # Create set of network IDs for links
    link_set = collect(keys(network.links))

    link_tot = Dict{}()
    link_length = Dict{}()
    # Dictionary to find node at other end of link from given node
    other_node = Dict{}()

    # Dictionary connecting nodes to zones
    for link in link_set
        link_length[link] = network.links[link].length
        dict_link_zones = filter(kv ->(kv[2].location == link), network.zones)
        link_zones = [(x.zone, x.fraction) for x in values(dict_link_zones)]
        dict_zones[link] = link_zones
    end

    # Populate the link related dictionaries
    for node in node_set
        # Links conected to a node (going to a node or emanating from a node)
        dict_links_tofrom_node = filter(kv -> ((kv[2].to == node) || (kv[2].from == node)), network.links)
        links_tofrom_node = sort(unique(keys(dict_links_tofrom_node)))
        link_tot[node] = links_tofrom_node
        for link in link_tot[node]
            if (network.links[link].from == node)
	            end_node = network.links[link].to
	        else
	            end_node = network.links[link].from
	        end
            # dictionary to find node at other end of link from given node
	        other_node[(node, link)] = end_node
	    end
    end

    maximum_number_of_links = maximum(length.(values(link_tot)))
    ################################# PATHWAYS ##################################################
    println("Processing Pathways")
    pathway_technology = Dict{Tuple{Pathway,Int64},Technology}()

    pathway_set = sort(unique(getfield.(keys(processlibrary.pathways), :pathway)))

    for (key,value) in processlibrary.pathways
	    pathway_technology[(key.pathway, key.stage)] = value.technology
    end

    pathway_stages = Dict{String, Int64}()
    pathway_unextended_storage_stages = Dict(p => Int64[] for p in pathway_set)
    pathway_extended_storage_stages = Dict(p => Int64[] for p in pathway_set)
    pathways_extended = Dict{String, Vector{Int64}}()
    pathways_unextended = Dict{String, Vector{Int64}}()
    pathway_stage_yield = Dict{Tuple{Pathway, Int64},Float64}()

    # Apportioning stages for each pathway into unextended, extended and other more detailed categories
    for pathway in pathway_set
	    pathway_keys = filter(k ->(k.pathway == pathway), keys(processlibrary.pathways))
        stage_set = sort(unique(getfield.(pathway_keys, :stage)))

        pathway_stages[pathway] = maximum(stage_set)

	    dict_pathway_extended = filter(kv ->(kv[1].pathway == pathway && kv[2].extended == true), processlibrary.pathways)
        dict_pathway_unextended = filter(kv ->(kv[1].pathway == pathway && kv[2].extended == false), processlibrary.pathways)

        # Extended pathway stages
        pathways_extended[pathway] = sort(unique(getfield.(keys(dict_pathway_extended), :stage)))
        # Unextended pathway stages
        pathways_unextended[pathway] = sort(unique(getfield.(keys(dict_pathway_unextended), :stage)))

        for stage in pathways_unextended[pathway]
            storage = storagelookup(processlibrary, pathway_technology[pathway, stage], total_year_set[1], 1.0, 0.0, "", "")
            if (storage.technology != "No Storage") && !isnothing(storage)
                push!(pathway_unextended_storage_stages[pathway], stage)
            end
        end

        for stage in pathways_extended[pathway]
            storage = storagelookup(processlibrary, pathway_technology[pathway, stage], total_year_set[1], 1.0, 0.0, "", "")
            if (storage.technology != "No Storage") && !isnothing(storage)
                push!(pathway_extended_storage_stages[pathway], stage)
            end
        end

        for stage in stage_set
            pathway_stage_yield[(pathway, stage)] = getfield(processlibrary.pathways[PathwayKey(pathway, stage)], :yield)
            pathway_stage_yield[(pathway, stage)] = 1.0
        end
    end

    pathways_unextended_max_length = maximum(maximum.(collect(values(pathways_unextended))))
    pathways_extended_max_length = maximum(maximum.(collect(values(pathways_extended))))

    technology_storage = Dict{Technology, Union{Storage,Nothing}}()

    for (key,value) in processlibrary.costs
	      technology_storage[(key.technology)] = value.storage
    end

    outputs_prices = joinpath("outputs", "hydrogen_prices_SERA_2.tsv")
    outputs_flow = joinpath("outputs", "flow_SERA_2.tsv")
    outputs_construction = joinpath("outputs", "construction_SERA_2.tsv")

    if isdir("outputs")
        if isfile(outputs_flow)
            rm(outputs_flow)
            rm(outputs_construction)
            rm(outputs_prices)
        end
    else
        mkdir("outputs")
    end

    rolling_production_available_cap = Dict{Tuple{String, String, String, String}, Float64}()
    rolling_unext_available_cap = Dict{Tuple{String, String, Int64, String, String}, Float64}()
    rolling_unext_storage_available_cap = Dict{Tuple{String, String, Int64, String, String}, Float64}()
    rolling_ext_available_cap = Dict{Tuple{String, String, Int64, String, String, String}, Float64}()
    rolling_ext_storage_available_cap = Dict{Tuple{String, String, Int64, String, String, String}, Float64}()
    iteration_ext_storage_cap = Dict{Tuple{String, String, Int64, String, String}, Float64}()
    previous_ext_storage_cap = Dict{Tuple{String, String, Int64, String, String}, Float64}()
    normalized_ext_storage_cap = Dict{Tuple{String, String, Int64, String, String}, Float64}()
    rolling_storage_level = Dict{Tuple{String, String, Int64, String}, Float64}()
    rolling_extended_storage_level = Dict{Tuple{String, String, String, Int64, String}, Float64}()

    infrastructure_id_set = String[]

    # Set initial flow for pathways
    init_flow = maximum(vcat(collect(values(demands_total)), collect(values(existings_capacity))))

    maximum_demand = maximum(values(demands_total))

    capacity_dict = Dict{String, Float64}()
    utilization_dict = Dict{String, Float64}()

    println("Creating Infrastructure options")
    @time begin
    for n in node_set
        for y in total_year_set
            if n in production_nodes
                for x in node_production_set[n]
                        infrastructureid = replace(string("INFR", '_', n, '_', x, '_', y), " " => "_")
                        push!(infrastructure_id_set, infrastructureid)
                        capacity, utilization = find_init_capacity_and_utilization(processlibrary, x, y)
                        capacity_dict[infrastructureid] = capacity * min_period_interval
                        utilization_dict["$(n)_$(x)_$(y)"] = utilization
                end
            end

            for p in pathway_set
                for s in pathways_unextended[p]
                        infrastructureid = replace(string("INFR", '_', n, '_', p, '_', s, '_', y), " " => "_")
                        push!(infrastructure_id_set, infrastructureid)
                        technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                        capacity, utilization = find_init_capacity_and_utilization(processlibrary, technology, y)
                        capacity_dict[infrastructureid] = capacity * min_period_interval
                        utilization_dict["$(n)_$(p)_$(s)_$(y)"] = utilization
                end

                for s in pathways_extended[p]
                        technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                        capacity, utilization = find_init_capacity_and_utilization(processlibrary, technology, y)
                        utilization_dict["$(n)_$(p)_$(s)_$(y)"] = utilization
                        for l in link_tot[n]
                            infrastructureid = replace(string("INFR", '_', n, '_', p,'_', s, '_', l, '_', y), " " => "_")
                            push!(infrastructure_id_set, infrastructureid)
                            technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                            capacity, utilization = find_init_capacity_and_utilization(processlibrary, technology, y)
                            capacity_dict[infrastructureid] = capacity * min_period_interval * 1.0
                        end
                end

                for s in pathway_unextended_storage_stages[p]
                        infrastructureid = replace(string("INFR", "_STORAGE_", n, '_', p, '_', s, '_', y), " " => "_")
                        push!(infrastructure_id_set, infrastructureid)
                        technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                        capacity, utilization = find_init_capacity_and_utilization(processlibrary, technology, y)
                        capacity_dict[infrastructureid] = capacity * min_period_interval
                        utilization_dict["$(n)_$(p)_$(s)_$(y)"] = utilization
                end
            end
        end
    end

    for n in node_set
        for y in total_year_set
            for t in total_year_period_set
                if n in production_nodes
                    for x in node_production_set[n]
                        rolling_production_available_cap[n, x, string(y), string(t[1])] = 0.0
                    end
                end

                for x in node_existings[n]
                    infrastructureid = replace(string("INFR", "_", n, '_', x), " " => "_")
                    push!(infrastructure_id_set, infrastructureid)
                    capacity_dict[infrastructureid] = existings_capacity[n, x, t]
                end

                for p in pathway_set
                    for s in pathways_unextended[p]
                        rolling_unext_available_cap[n, p, s, string(y), string(t[1])] = 0.0
                    end

                    for s in pathways_extended[p]
                            for l in link_tot[n]
                                rolling_ext_available_cap[n, p, s, l, string(y), string(t[1])] = 0.0
                                if s in pathway_extended_storage_stages[p]
                                    rolling_ext_storage_available_cap[n, p, s, l, string(y), string(t[1])] = 0.0
                                    previous_ext_storage_cap[n, p, s, l, string(y)] = 0.0
                                    normalized_ext_storage_cap[n, p, s, l, string(y)] = 1.0
                                    rolling_extended_storage_level[n, l, p, s, string(y)] = 0.0
                                end
                            end
                    end

                    for s in pathway_unextended_storage_stages[p]
                        rolling_unext_storage_available_cap[n, p, s, string(y), string(t[1])] = 0.0
                        rolling_storage_level[n, p, s, string(y)] = 0.0
                    end
                end
            end
        end
    end

    end

    production_data = Dict{Tuple{Int64, String, String, String}, Construction}()
    production_annualization = Dict{Tuple{Int64, String, String, String}, Float64}()

    pathway_unext_flow_data = Dict{Tuple{Int64, String, String, Int64, String}, Construction}()
    pathway_unext_annualization = Dict{Tuple{Int64, String, String, Int64, String}, Float64}()

    pathway_ext_link_data = Dict{Tuple{Int64, String, String, Int64, String, String}, Construction}()
    pathway_ext_annualization = Dict{Tuple{Int64, String, String, Int64, String, String}, Float64}()

    pathway_unext_storage_data = Dict{Tuple{Int64, String, String, Int64, String}, Construction}()
    pathway_storage_annualization = Dict{Tuple{Int64, String, String, Int64, String}, Float64}()

    production_price_per_kg_of_hydrogen = Dict{Tuple{String, String, Tuple{Year,PeriodName}}, Float64}()
    existings_price_per_kg_of_hydrogen = Dict{Tuple{String, String, Tuple{Year,PeriodName}}, Float64}()
    pathway_unextended_price_per_kg_of_hydrogen = Dict{Tuple{String, String, Int64, Tuple{Year,PeriodName}}, Float64}()
    pathway_extended_price_per_kg_of_hydrogen = Dict{Tuple{String, String, Int64, String, Tuple{Year,PeriodName}}, Float64}()
    storage_price_per_kg_of_hydrogen = Dict{Tuple{String, String, Int64, Tuple{Year,PeriodName}}, Float64}()
    available_extended_storage_cap = Dict{Tuple{String, String, String, Int64, Int64, Int64}, Float64}()

    existings_data = Dict{Tuple{String, String, Tuple{Year,PeriodName}}, Existing}()

    fixed_cost_dict = Dict{String, Float64}()

    for current_year in start_year:rolling_step:total_year_set[end] - rolling_window + rolling_step

        start_index = findfirst(x -> x == current_year, total_year_set)
        if !(isnothing(start_index))
            end_index = min(start_index + rolling_window - 1, min(length(total_year_set), findfirst(x -> x == last_year, total_year_set)))
            if end_index >= start_index'
                printstyled("\nOPTIMIZING FROM $(current_year) TO $(total_year_set[end_index]) \n\n"; color = :yellow)
                year_set = total_year_set[start_index:end_index]
                year_period_set = filter(x -> x[1] in year_set, total_year_period_set)
                year_set_string = string.(year_set)

                model = Model(Xpress_optimizer)

                ############################################# CREATE OPTIMIZATION MODEL ########################################################################

                println("STARTING OPTIMIZATION MODEL CREATION")
                ##### VARIABLES #################

                println("CREATING VARIABLES AND EXPRESSIONS")
                @time begin
                # Planning Variables:
                if build_new_production
                    @variable(model, new_production_cap[n in production_nodes, x in node_production_set[n], t in year_set] >= 0.0)
                else
                    @variable(model, new_production_cap[n in production_nodes, x in node_production_set[n], t in year_set] == 0.0)
                end
                @variable(model, pathway_unextended_cap[n in node_set, p in pathway_set, s in pathways_unextended[p], t in year_set] >= 0.0)
                @variable(model, pathway_extended_cap[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], t in year_set] >= 0.0)

                if storage_available
                    @variable(model, pathway_storage_cap[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], t in year_set] >= 0.0)

                    # Variables for hydrogen flowing into storage in extended stages
                    @variable(model, extended_stage_stor_in[n in node_set, link in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in total_year_set, t in year_period_set] >= 0.0)

                    # Variables for hydrogen flowing into storage in extended stages
                    @variable(model, extended_stage_stor_out[n in node_set, link in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in total_year_set, t in year_period_set] >= 0.0)

                else
                    @variable(model, pathway_storage_cap[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], t in year_set] == 0.0)

                    # Variables for hydrogen flowing into storage in extended stages
                    @variable(model, extended_stage_stor_in[n in node_set, link in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in total_year_set, t in year_period_set] == 0.0)

                    # Variables for hydrogen flowing into storage in extended stages
                    @variable(model, extended_stage_stor_out[n in node_set, link in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in total_year_set, t in year_period_set] == 0.0)

                end

                # Operation Variables:
                # Production variable
                @variable(model, new_production[n in production_nodes, x in node_production_set[n], y in total_year_set, t in year_period_set] >= 0.0)

                # Existing variable
                @variable(model, existing_production[n in node_set, x in node_existings[n], t in year_period_set] >= 0.0)

                # Flow variables for all pathway stages
                @variable(model, pathway_stage_flow[n in node_set, p in pathway_set, s in 1:pathway_stages[p], y in total_year_set, t in year_period_set] >= 0.0)

                # Variables for hydrogen flowing into storage at unextended stages
                @variable(model, pathway_stage_stor_in[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], y in total_year_set, t in year_period_set] >= 0.0)

                # Variables for hydrogen flowing out of storage at unextended stages
                @variable(model, pathway_stage_stor_out[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], y in total_year_set, t in year_period_set] >= 0.0)

                # Variables for hydrogen stored at unextended stages
                @variable(model, pathway_stage_stor_level[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], y in total_year_set, t in year_period_set] >= 0.0)

                # Variables for hydrogen flowing from node n across link l
                @variable(model, link_flow[n in node_set, link in link_tot[n], p in pathway_set, s in pathways_extended[p], y in total_year_set, t in year_period_set] >= 0.0)

                # Dummy Variables
                @variable(model, empty == 0.0)

                ##### EXPRESSIONS #################

                # Expression that consumption equals total demands
                @expression(model, consumption[n in demands_node_set, t in year_period_set], demands_total[(n, t)])

                # Aggregate production expression
                @expression(model, aggregate_production[n in node_set, t in year_period_set], 1 * empty)

                # Expression for net hydrogen flow at each unextended stage of the pathway
                @expression(model, stage_net_flow[n in node_set, p in pathway_set, s in pathways_unextended[p], y in total_year_set, t in year_period_set], pathway_stage_yield[p, s] * pathway_stage_flow[n, p, s, y, t])

                # Expression for total available capacity based on new production added each year (where y is the construction year)
                JuMP.@expression(model, available_production_cap[n in production_nodes, x in node_production_set[n], y in total_year_set, t in year_set], rolling_production_available_cap[n, x, string(y), string(t)] * utilization_dict["$(n)_$(x)_$(y)"] + empty)

                # Expression for total available capacity based on new unextended pathways added each year (where y is the construction year)
                JuMP.@expression(model, available_unextended_cap[n in node_set, p in pathway_set, s in pathways_unextended[p], y in total_year_set, t in year_set], rolling_unext_available_cap[n, p, s, string(y), string(t)] * utilization_dict["$(n)_$(p)_$(s)_$(y)"] + empty)

                # Expression for total available capacity based on new storage pathways added each year (where y is the construction year)
                JuMP.@expression(model, available_storage_cap[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], y in total_year_set, t in year_set], rolling_unext_storage_available_cap[n, p, s, string(y), string(t)] * utilization_dict["$(n)_$(p)_$(s)_$(y)"] + empty)

                # Expression for total available capacity based on new extended pathways added each year (where y is the construction year)
                JuMP.@expression(model, available_extended_cap[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in total_year_set, t in year_set], rolling_ext_available_cap[n, p, s, l, string(y), string(t)] * utilization_dict["$(n)_$(p)_$(s)_$(y)"] + empty)

                # Expression for hydrogen stored in extended stages

                JuMP.@expression(model, extended_stage_stor_level[n in node_set, link in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in total_year_set, t in year_period_set], 1.0 * empty)
                for n in node_set
                    for l in link_tot[n]
                        for p in pathway_set
                            for s in pathway_extended_storage_stages[p]
                                for y in total_year_set
                                    for t in year_set
                                        available_extended_storage_cap[n, l, p, s, y, t] = rolling_ext_storage_available_cap[n, p, s, l, string(y), string(t)]
                                    end
                                end
                            end
                        end
                    end
                end

                end
                println("POPULATING EXPRESSIONS")

                @time begin
                for n in node_set
                    for t in year_period_set
                        for x in node_existings[n]
                            add_to_expression!(aggregate_production[n, t], existing_production[n, x, t])
                        end
                        if (n in production_nodes) && build_new_production
                            add_to_expression!(aggregate_production[n, t], sum(new_production[n, x, y, t] for x in node_production_set[n], y in filter(x -> x <= t[1], total_year_set)))
                        end
                    end

                    for p in pathway_set
                        for s in pathway_unextended_storage_stages[p]
                            for y in total_year_set
                                for t in year_period_set
                                    add_to_expression!(stage_net_flow[n, p, s, y, t], pathway_stage_stor_out[n, p, s, y, t] - pathway_stage_stor_in[n, p, s, y, t])
                                end
                            end
                        end
                    end

                    for t in year_set
                        for y in current_and_prev_years(t, year_set)
                            if n in production_nodes
                                for x in node_production_set[n]
                                    if y + processcosts(processlibrary, x, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                                        JuMP.add_to_expression!(available_production_cap[n, x, y, t], new_production_cap[n, x, y] * utilization_dict["$(n)_$(x)_$(y)"])
                                    end
                                end
                            end

                            for p in pathway_set
                                for s in pathways_unextended[p]
                                    technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                    if y + processcosts(processlibrary, technology, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                                        JuMP.add_to_expression!(available_unextended_cap[n, p, s, y, t], pathway_unextended_cap[n, p, s, y] * utilization_dict["$(n)_$(p)_$(s)_$(y)"])
                                    end
                                end

                                for s in pathway_unextended_storage_stages[p]
                                    technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                    storage = storagelookup(processlibrary, technology, y, 0.0, 0.0, "", n)
                                    if y + processcosts(processlibrary, storage.technology, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                                        JuMP.add_to_expression!(available_storage_cap[n, p, s, y, t], pathway_storage_cap[n, p, s, y] * utilization_dict["$(n)_$(p)_$(s)_$(y)"])
                                    end
                                end

                                for s in pathways_extended[p]
                                    for l in link_tot[n]
                                        technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                        if y + processcosts(processlibrary, technology, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                                            JuMP.add_to_expression!(available_extended_cap[n, l, p, s, y, t], pathway_extended_cap[n, l, p, s, y] * utilization_dict["$(n)_$(p)_$(s)_$(y)"])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
                end

                 # Expression for total available storage capacity.
                 JuMP.@expression(model, total_storage_capacity[t in year_set],
                                        sum(available_storage_cap[n, p, s, y, t] for n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], y in total_year_set) +
                                        sum(rolling_ext_storage_available_cap[n, p, s, l, string(y), string(t)] for n in node_set, l in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in total_year_set))

                final_period = last(period_set)
                maximum_capacity_technologies = keys(maximum_capacity_constraints)

                println("CREATING CONSTRAINTS")
                @time begin
                for n in node_set
                    for (t_idx, t) in enumerate(year_period_set)

                        @constraint(model, aggregate_production[n, t] >= 0.0)
                        @constraint(model, aggregate_production[n, t] == sum(pathway_stage_flow[n, p, 1, y, t]
                                                                        for p in pathway_set, y in total_year_set))

                        if n in production_nodes
                            for x in node_production_set[n]
                                utilization_requirement = 0.0
                                if x in keys(utilization_constraints)
                                    utilization_requirement = utilization_constraints[x]
                                end

                                for y in total_year_set
                                    @constraint(model, new_production[n, x, y, t] <= available_production_cap[n, x, y, t[1]] * periods[t[2]])
                                    if t[2] == final_period
                                        @constraint(model, sum(new_production[n, x, y, (t[1], h)] for h in period_set) >= utilization_requirement * available_production_cap[n, x, y, t[1]])
                                    end
                                    if x in maximum_capacity_technologies
                                            @constraint(model, new_production[n, x, y, t] <= maximum_capacity_constraints[x] * periods[t[2]])
                                    end
                                end
                            end

                            if n in onsite_nodes
                                @constraint(model, sum(new_production[n, x, y, t] for x in node_production_set[n], y in total_year_set) <= demands_total[n, t])
                            end
                        end

                        for x in node_existings[n]
                            @constraint(model, existing_production[n, x, t] <= existings_capacity[(n, x, t)])

                            utilization_requirement = 0.0
                            if x in keys(utilization_constraints)
                                utilization_requirement = utilization_constraints[x]
                            end
                            if t[2] == final_period
                                @constraint(model, sum(existing_production[n, x, (t[1], h)] for h in period_set) >= utilization_requirement * sum(existings_capacity[(n, x, (t[1], h))] for h in period_set))
                            end
                        end

                        for p in pathway_set
                            for s in pathways_unextended[p]
                                technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                utilization_requirement = 0.0
                                if technology in keys(utilization_constraints)
                                    utilization_requirement = utilization_constraints[technology]
                                end

                                if s != last(pathways_unextended[p])
                                    @constraint(model, sum(stage_net_flow[n, p, s, y, t] for y in total_year_set) == sum(pathway_stage_flow[n, p, s + 1, y, t] for y in total_year_set))
                                end

                                for y in total_year_set

                                    @constraint(model, pathway_stage_flow[n, p, s, y, t] <= available_unextended_cap[n, p, s, y, t[1]] * periods[t[2]])
                                    if t[2] == final_period
                                        @constraint(model, sum(pathway_stage_flow[n, p, s, y, (t[1], h)] for h in period_set) >= utilization_requirement * available_unextended_cap[n, p, s, y, t[1]])
                                    end

                                    if technology in maximum_capacity_technologies
                                        @constraint(model, pathway_stage_flow[n, p, s, y, t] <= maximum_capacity_constraints[technology] * periods[t[2]])
                                    end
                                end
                            end

                            for s in pathway_unextended_storage_stages[p]
                                technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                utilization_requirement = 0.0
                                if technology in keys(utilization_constraints)
                                    utilization_requirement = utilization_constraints[technology]
                                end
                                for y in total_year_set
                                    if t_idx == 1
                                        @constraint(model, pathway_stage_stor_level[n, p, s, y, t] == rolling_storage_level[n, p, s, string(y)] + pathway_stage_stor_in[n, p, s, y, t] - pathway_stage_stor_out[n, p, s, y, t])
                                    else
                                        @constraint(model, pathway_stage_stor_level[n, p, s, y, t] == pathway_stage_stor_level[n, p, s, y, prev(t, period_set, year_set)] + pathway_stage_stor_in[n, p, s, y, t] - pathway_stage_stor_out[n, p, s, y, t])
                                    end

                                    @constraint(model, pathway_stage_stor_level[n, p, s, y, t] <= available_storage_cap[n, p, s, y, t[1]])

                                end
                            end

                            for s in pathways_extended[p]
                                technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                utilization_requirement = 0.0
                                if technology in keys(utilization_constraints)
                                    utilization_requirement = utilization_constraints[technology]
                                end
                                if s in pathway_extended_storage_stages[p]

                                    @constraint(model,  sum(pathway_stage_flow[n, p, s, y, t] +
                                                            sum((pathway_stage_yield[p, s] * link_flow[other_node[(n, l)], l, p, s, y, t]) -
                                                            link_flow[n, l, p, s, y, t] +
                                                            extended_stage_stor_out[other_node[(n, l)], l, p, s, y, t] -
                                                            extended_stage_stor_in[other_node[(n, l)], l, p, s, y, t]
                                                              for l in link_tot[n]) for y in total_year_set)  ==
                                                            sum(pathway_stage_flow[n, p, s + 1, y, t] for y in total_year_set) )

                                    for y in total_year_set
                                        for l in link_tot[n]

                                            if t_idx == 1
                                                JuMP.add_to_expression!(extended_stage_stor_level[n, l, p, s, y, t], rolling_extended_storage_level[n, l, p, s, string(y)] + extended_stage_stor_in[n, l, p, s, y, t] - extended_stage_stor_out[n, l, p, s, y, t])
                                            else
                                                JuMP.add_to_expression!(extended_stage_stor_level[n, l, p, s, y, t], extended_stage_stor_level[n, l, p, s, y, prev(t, period_set, year_set)] + extended_stage_stor_in[n, l, p, s, y, t] - extended_stage_stor_out[n, l, p, s, y, t])
                                            end

                                            @constraint(model, extended_stage_stor_level[n, l, p, s, y, t] >= 0.0)

                                            @constraint(model, extended_stage_stor_in[n, l, p, s, y, t] <= link_flow[n, l, p, s, y, t])

                                            @constraint(model, link_flow[n, l, p, s, y, t] <= available_extended_cap[n, l, p, s, y, t[1]] * periods[t[2]])
                                            if t[2] == final_period
                                                @constraint(model, sum(link_flow[n, l, p, s, y, (t[1], h)] for h in period_set) >= utilization_requirement * available_extended_cap[n, l, p, s, y, t[1]])
                                            end
                                            if technology in maximum_capacity_technologies
                                                @constraint(model, link_flow[n, l, p, s, y, t] <= maximum_capacity_constraints[technology] * periods[t[2]])
                                            end
                                        end
                                    end

                                else
                                    @constraint(model,  sum(pathway_stage_flow[n, p, s, y, t] +
                                                    sum((pathway_stage_yield[p, s] * link_flow[other_node[(n, l)], l, p, s, y, t]) -
                                                    link_flow[n, l, p, s, y, t]
                                                        for l in link_tot[n]) for y in total_year_set)  ==
                                                    sum(pathway_stage_flow[n, p, s + 1, y, t] for y in total_year_set) )

                                    for y in total_year_set
                                        for l in link_tot[n]

                                            @constraint(model, link_flow[n, l, p, s, y, t] <= available_extended_cap[n, l, p, s, y, t[1]] * periods[t[2]])
                                            if t[2] == final_period
                                                @constraint(model, sum(link_flow[n, l, p, s, y, (t[1], h)] for h in period_set) >= utilization_requirement * available_extended_cap[n, l, p, s, y, t[1]])
                                            end
                                            if technology in maximum_capacity_technologies
                                                @constraint(model, link_flow[n, l, p, s, y, t] <= maximum_capacity_constraints[technology] * periods[t[2]])
                                            end
                                        end
                                    end
                                end

                            end
                        end

                    end
                end

                @constraint(model, extended_storage_level_limit[n in node_set, l in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in filter(j -> j in year_set, total_year_set), t in year_period_set],
                            extended_stage_stor_level[n, l, p, s, y, t] <= (rolling_ext_storage_available_cap[n, p, s, l, string(y), string(t[1])] + pathway_extended_cap[n, l, p, s, y] * normalized_ext_storage_cap[n, p, s, l, string(y)]) * storage_available)

                @constraint(model, extended_storage_level_limit_existing[n in node_set, l in link_tot[n], p in pathway_set, s in pathway_extended_storage_stages[p], y in filter(j -> !(j in year_set), total_year_set), t in year_period_set],
                        extended_stage_stor_level[n, l, p, s, y, t] <= rolling_ext_storage_available_cap[n, p, s, l, string(y), string(t[1])] * storage_available)

                # Sum of net H2 flowing through the last stage of all pathways at node n should be equal to the total H2 demand at node n:
                @constraint(model, last_stage_to_consumption[n in demands_node_set, t in year_period_set],
                sum(stage_net_flow[n, p, pathway_stages[p], y, t] for p in pathway_set, y in total_year_set) >= consumption[n, t])

                # minimum percentage of H2 demand should be met by eligible technologies

                if !isempty(clean_production_constraints)
                    for constr in clean_production_constraints
                        if constr.year in year_set
                            @constraint(model,
                            sum(new_production[n, x, y, (constr.year, h)]
                                for n in production_nodes, x in filter(j -> j in constr.eligible_technologies, node_production_set[n]), y in total_year_set, h in period_set) +
                            sum(existing_production[n, x, (constr.year, h)]
                                for n in existings_node_set, x in filter(y -> y in constr.eligible_technologies, node_existings[n]), h in period_set) >=
                            sum(consumption[n, (constr.year, h)] for n in demands_node_set, h in period_set) * constr.requirement)
                        end
                    end
                end

                # minimum production capacity requirement to be met by eligible technologies
                if !isempty(production_capacity_constraints)
                    if production_capacity_switch == "Nodal"
                        for constr in production_capacity_constraints
                            if constr.year in year_set && constr.location in node_set
                                @constraint(model,
                                sum(available_production_cap[constr.location, x, y, constr.year]
                                    for x in filter(y -> y in constr.eligible_technologies, node_production_set[constr.location]), y in total_year_set) +
                                sum(existing_capacity[constr.location, x, constr.year]
                                    for x in filter(y -> y in constr.eligible_technologies, existings_set)) >=
                                constr.capacity)
                            end
                        end
                    elseif production_capacity_switch == "SystemWide"
                        systemwide_cap_constraint = Dict{Tuple{Year, Vector{Technology}}, Float64}()
                        for constr in production_capacity_constraints
                            systemwide_cap_constraint[(constr.year, sort(constr.eligible_technologies))] = 0.0
                        end
                        for constr in production_capacity_constraints
                            systemwide_cap_constraint[(constr.year, sort(constr.eligible_technologies))] += constr.capacity
                        end
                        for constr in keys(systemwide_cap_constraint)
                            if constr[1] in year_set
                                println(constr)
                                println(systemwide_cap_constraint[constr])
                                @constraint(model,
                                sum(available_production_cap[n, x, y, constr[1]] 
                                    for n in production_nodes, x in filter(w -> w in constr[2], node_production_set[n]), y in total_year_set) + 
                                sum(existing_capacity[n, x, constr[1]]
                                    for n in existings_node_set, x in filter(w -> w in constr[2], node_existings[n])) >= 
                                systemwide_cap_constraint[constr])
                             end
                        end
                    end
                end

                # minimum storage capacity requirement
                if !isempty(storage_capacity_constraints)
                    storage_constr_years = filter(x -> x in year_set, keys(storage_capacity_constraints))
                    @constraint(model, total_storage_capacity_constr[y in storage_constr_years], total_storage_capacity[y] >= storage_capacity_constraints[y])
                end

                end

                println("POPULATED CORE OPTIMIZATION MODEL, STARTING ITERATIONS")

                # Objective Function of the model requires cost data which needs to be updated iteratively.
                #Start Iterations
                i = 1
                ϵ = 1e8
                prev_objective = 0.0

                while i <= niter && ϵ > ϵ_max
                    println("READING COST DATA FOR ITERATION $(i)")
                    @time begin
                    for n in node_set
                        for t in year_set
                            remaining_years = filter(y -> y >= t, year_set)
                            if n in production_nodes
                                for x in node_production_set[n]
                                    infrastructureid = replace(string("INFR", '_', n, '_', x, '_', t), " " => "_")
                                    #capacity = flowDF[findfirst(x -> x == infrastructureid, flowDF[!, :Infrastructure_ID]), :Flow]
                                    capacity = capacity_dict[infrastructureid]
                                    distance = 0.0
                                    production_data[i, n, x, string(t)] = processcosts(processlibrary, x, t, capacity, distance, infrastructureid, n)
                                    production_inputs = processinputs(processlibrary, x, t, capacity, distance)
                                    production_outputs = processoutputs(processlibrary, x, t, capacity, distance)
                                    input_materials = keys(production_inputs)
                                    for h in period_set
                                        production_price_per_kg_of_hydrogen[n, x, (t, h)] = 0.0
                                        for material in input_materials
                                            # For this particular material, look up consumption/(kg of H2)
                                            quantity = production_inputs[material]
                                            price_material = 0.0
                                            for z in 1:length(dict_zones[n])
                                                zone = dict_zones[n][z][1]
                                                fraction = dict_zones[n][z][2]
                                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                            end
                                            production_price_per_kg_of_hydrogen[n, x, (t, h)] += price_material * quantity
                                        end

                                        output_materials = keys(production_outputs)
                                        for material in output_materials
                                            # For this particular material, look up output/(kg of H2)
                                            quantity = production_outputs[material]
                                            price_material = 0.0
                                            for z in 1:length(dict_zones[n])
                                                zone = dict_zones[n][z][1]
                                                fraction = dict_zones[n][z][2]
                                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                            end
                                            production_price_per_kg_of_hydrogen[n, x, (t, h)] += price_material * quantity
                                        end
                                        if isinf(production_price_per_kg_of_hydrogen[n, x, (t, h)]) || production_price_per_kg_of_hydrogen[n, x, (t, h)] > 1e10
                                            @constraint(model, new_production_cap[n, x, t] == 0.0)
                                            production_price_per_kg_of_hydrogen[n, x, (t, h)] = 0.0
                                        end

                                    end

                                    if iszero(discount_rate)
                                        production_annualization[i, n, x, string(t)] = 1 / production_data[i, n, x, string(t)].lifetime
                                    else
                                        production_annualization[i, n, x, string(t)] = discount_rate / (1 - (1 + discount_rate) ^ (-(production_data[i, n, x, string(t)].lifetime)))
                                    end

                                end

                            end


                            for p in pathway_set
                                for s in pathways_unextended[p]
                                    infrastructureid = replace(string("INFR", '_', n, '_', p, '_', s, '_', t), " " => "_")
                                    capacity = capacity_dict[infrastructureid] / min_period_interval
                                    distance = 0.0
                                    technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                    pathway_unext_flow_data[i, n, p, s, string(t)] = processcosts(processlibrary, technology, t, capacity, distance, infrastructureid, n)

                                    pathway_unext_flow_inputs = processinputs(processlibrary, technology, t, capacity, distance)
                                    pathway_unext_flow_outputs = processoutputs(processlibrary, technology, t, capacity, distance)

                                    input_materials = keys(pathway_unext_flow_inputs)
                                    for h in period_set
                                        pathway_unextended_price_per_kg_of_hydrogen[n, p, s, (t, h)] = 0.0
                                        for material in input_materials
                                            # For this particular material, look up consumption/(kg of H2)
                                            quantity = pathway_unext_flow_inputs[material]
                                            price_material = 0.0
                                            for z in 1:length(dict_zones[n])
                                                zone = dict_zones[n][z][1]
                                                fraction = dict_zones[n][z][2]
                                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                            end
                                            pathway_unextended_price_per_kg_of_hydrogen[n, p, s, (t, h)] += price_material * quantity
                                        end

                                        output_materials = keys(pathway_unext_flow_outputs)
                                        for material in output_materials
                                            # For this particular material, look up output/(kg of H2)
                                            quantity = pathway_unext_flow_outputs[material]
                                            price_material = 0.0
                                            for z in 1:length(dict_zones[n])
                                                zone = dict_zones[n][z][1]
                                                fraction = dict_zones[n][z][2]
                                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                            end
                                            pathway_unextended_price_per_kg_of_hydrogen[n, p, s, (t, h)] += price_material * quantity
                                        end

                                        if isinf(pathway_unextended_price_per_kg_of_hydrogen[n, p, s, (t, h)]) || pathway_unextended_price_per_kg_of_hydrogen[n, p, s, (t, h)] > 1e10
                                            @constraint(model, pathway_unextended_cap[n, p, s, t] == 0.0)
                                            pathway_unextended_price_per_kg_of_hydrogen[n, p, s, (t, h)] = 0.0
                                        end
                                    end
                                    if iszero(discount_rate)
                                        pathway_unext_annualization[i, n, p, s, string(t)] = 1 / pathway_unext_flow_data[i, n, p, s, string(t)].lifetime
                                    else
                                        pathway_unext_annualization[i, n, p, s, string(t)] = discount_rate / (1 - (1 + discount_rate) ^ (-(pathway_unext_flow_data[i, n, p, s, string(t)].lifetime)))
                                    end

                                    storage = storagelookup(processlibrary, technology, t, capacity, distance, infrastructureid, n)
                                    if ((storage.technology != "No Storage") && !isnothing(storage))
                                        infrastructureid = replace(string("INFR", "_STORAGE_", n, '_', p, '_', s, '_', t), " " => "_")
                                        capacity = capacity_dict[infrastructureid] / min_period_interval
                                        pathway_unext_storage_data[i, n, p, s, string(t)] = processcosts(processlibrary, storage.technology, t, capacity, distance, infrastructureid, n)

                                        #@constraint(model, pathway_storage_cap[n, p, s, t] <= pathway_unext_storage_data[i, n, p, s, string(t)].storagecapacity * storage_available)

                                        pathway_storage_inputs = processinputs(processlibrary, storage.technology, t, capacity, distance)
                                        pathway_storage_outputs = processoutputs(processlibrary, storage.technology, t, capacity, distance)

                                        input_materials = keys(pathway_storage_inputs)
                                        for h in period_set
                                            storage_price_per_kg_of_hydrogen[n, p, s, (t, h)] = 0.0
                                            for material in input_materials
                                            # For this particular material, look up consumption/(kg of H2)
                                                quantity = pathway_storage_inputs[material]
                                                price_material = 0.0
                                                for z in 1:length(dict_zones[n])
                                                    zone = dict_zones[n][z][1]
                                                    fraction = dict_zones[n][z][2]
                                                    price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                                end
                                                storage_price_per_kg_of_hydrogen[n, p, s, (t, h)] += price_material * quantity
                                            end
                                            output_materials = keys(pathway_storage_outputs)
                                            for material in output_materials
                                            # For this particular material, look up production/(kg of H2)
                                                quantity = pathway_storage_outputs[material]
                                                price_material = 0.0
                                                for z in 1:length(dict_zones[n])
                                                    zone = dict_zones[n][z][1]
                                                    fraction = dict_zones[n][z][2]
                                                    price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                                end
                                                storage_price_per_kg_of_hydrogen[n, p, s, (t, h)] += price_material * quantity
                                            end

                                            if isinf(storage_price_per_kg_of_hydrogen[n, p, s, (t, h)]) || storage_price_per_kg_of_hydrogen[n, p, s, (t, h)] > 1e10
                                                @constraint(model, pathway_storage_cap[n, p, s, t] == 0.0)
                                                storage_price_per_kg_of_hydrogen[n, p, s, (t, h)] = 0.0
                                            end
                                        end

                                        if iszero(discount_rate)
                                            pathway_storage_annualization[i, n, p, s, string(t)] = 1 / pathway_unext_storage_data[i, n, p, s, string(t)].lifetime
                                        else
                                            pathway_storage_annualization[i, n, p, s, string(t)] = discount_rate / (1 - (1 + discount_rate) ^ (-(pathway_unext_storage_data[i, n, p, s, string(t)].lifetime)))
                                        end
                                    end

                                end

                                for s in pathways_extended[p]
                                for l in link_tot[n]
                                    infrastructureid = replace(string("INFR", '_', n, '_', p,'_', s, '_', l, '_', t), " " => "_")
                                    capacity = capacity_dict[infrastructureid] / min_period_interval
                                    distance = link_length[l]

                                    technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)
                                    pathway_ext_link_data[i, n, p, s, l, string(t)] = processcosts(processlibrary, technology, t, capacity, distance, infrastructureid, l)
                                    pathway_ext_link_inputs = processinputs(processlibrary, technology, t, capacity, distance)
                                    pathway_ext_link_outputs = processoutputs(processlibrary, technology, t, capacity, distance)

                                    normalized_ext_storage_cap[n, p, s, l, string(t)] = 0.0

                                    if s in pathway_extended_storage_stages[p]
                                        normalized_ext_storage_cap[n, p, s, l, string(t)] = pathway_ext_link_data[i, n, p, s, l, string(t)].storagecapacity /  pathway_ext_link_data[i, n, p, s, l, string(t)].nameplate
                                        iteration_ext_storage_cap[n, p, s, l, string(t)] = normalized_ext_storage_cap[n, p, s, l, string(t)] * capacity
                                        for y in remaining_years
                                            available_extended_storage_cap[n, l, p, s, t, y] += iteration_ext_storage_cap[n, p, s, l, string(t)] - previous_ext_storage_cap[n, p, s, l, string(t)]
                                            if available_extended_storage_cap[n, l, p, s, t, y] < 1e-2
                                                available_extended_storage_cap[n, l, p, s, t, y] = 0.0
                                            end
                                        end

                                        previous_ext_storage_cap[n, p, s, l, string(t)] = deepcopy(iteration_ext_storage_cap[n, p, s, l, string(t)])

                                    end

                                    if pathway_ext_link_data[i, n, p, s, l, string(t)].storagecapacity < 0
                                        #@info "Storage capacity for linepack of link $(l) is less than 0. Will not build any linepack pipeline."
                                        @constraint(model, pathway_extended_cap[n, l, p, s, t] == 0.0)
                                        for y in remaining_years
                                            available_extended_storage_cap[n, l, p, s, t, y] = 0.0
                                        end
                                        pathway_ext_link_data[i, n, p, s, l, string(t)].capitalcost = 0.0
                                        pathway_ext_link_data[i, n, p, s, l, string(t)].fixedcost = 0.0
                                    end

                                    if pathway_ext_link_data[i, n, p, s, l, string(t)].capitalcost < 0
                                        #@info "Capital cost for linepack of link $(l) is less than 0. Will not build any linepack pipeline."
                                        @constraint(model, pathway_extended_cap[n, l, p, s, t] == 0.0)
                                        for y in remaining_years
                                            available_extended_storage_cap[n, l, p, s, t, y] = 0.0
                                        end
                                        pathway_ext_link_data[i, n, p, s, l, string(t)].capitalcost = 0.0
                                        pathway_ext_link_data[i, n, p, s, l, string(t)].fixedcost = 0.0
                                    end

                                    if (occursin("pipeline", lowercase(technology))) && (distance < min_pipeline_distance)
                                        @constraint(model, pathway_extended_cap[n, l, p, s, t] == 0.0)
                                       # @info "Length of link $(l) is less than $(min_pipeline_distance). Will not build any $(technology) on this link."
                                    end

                                    input_materials = keys(pathway_ext_link_inputs)
                                    for h in period_set

                                        pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, (t, h)] = 0.0
                                        for material in input_materials
                                            # For this particular material, look up consumption/(kg of H2)
                                            quantity = pathway_ext_link_inputs[material]
                                            price_material = 0.0
                                            for z in 1:length(dict_zones[l])
                                                zone = dict_zones[l][z][1]
                                                fraction = dict_zones[l][z][2]
                                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                            end
                                            pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, (t, h)] += price_material * quantity
                                        end

                                        output_materials = keys(pathway_ext_link_outputs)
                                        for material in output_materials
                                            # For this particular material, look up production/(kg of H2)
                                            quantity = pathway_ext_link_outputs[material]
                                            price_material = 0.0
                                            for z in 1:length(dict_zones[l])
                                                zone = dict_zones[l][z][1]
                                                fraction = dict_zones[l][z][2]
                                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, h, quantity ,true, true)
                                            end
                                            pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, (t, h)] += price_material * quantity
                                        end

                                        if isinf(pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, (t, h)]) || pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, (t, h)] > 1e10
                                            @constraint(model, pathway_extended_cap[n, l, p, s, t] == 0.0)
                                            pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, (t, h)] = 0.0
                                        end
                                    end

                                    if iszero(discount_rate)
                                        pathway_ext_annualization[i, n, p, s, l, string(t)] = 1 / pathway_ext_link_data[i, n, p, s, l, string(t)].lifetime
                                    else
                                        pathway_ext_annualization[i, n, p, s, l, string(t)] = discount_rate / (1 - (1 + discount_rate) ^ (-(pathway_ext_link_data[i, n, p, s, l, string(t)].lifetime)))
                                    end

                                end
                                end

                            end

                        end

                        for t in total_year_period_set
                            for x in node_existings[n]
                                capacity = existings_capacity[(n, x, t)]
                                cost = existings_cost[n, x, t]
                                location = n
                                technology = x
                                year = t[1]
                                period = t[2]
                                distance = 0.0
                                yield = existings_yield[n, x, t]
                                existings_data[n, x, t] = Existing(location, technology, year, period, capacity, yield, cost)

                                existings_inputs = processinputs(processlibrary, x, year, capacity, distance)
                                existings_outputs = processoutputs(processlibrary, x, year, capacity, distance)
                                input_materials = keys(existings_inputs)
                                existings_price_per_kg_of_hydrogen[n, x, t] = 0.0
                                for material in input_materials
                                    # For this particular material, look up consumption/(kg of H2)
                                    quantity = existings_inputs[material]
                                    price_material = 0.0
                                    for z in 1:length(dict_zones[n])
                                        zone = dict_zones[n][z][1]
                                        fraction = dict_zones[n][z][2]
                                        price_material = price_material + fraction * pricelookup(prices, usage, material, zone, year, period, quantity ,true, true)
                                    end
                                    existings_price_per_kg_of_hydrogen[n, x, t] += price_material * quantity
                                end

                                output_materials = keys(existings_outputs)
                                for material in output_materials
                                    # For this particular material, look up production/(kg of H2)
                                    quantity = existings_outputs[material]
                                    price_material = 0.0
                                    for z in 1:length(dict_zones[n])
                                        zone = dict_zones[n][z][1]
                                        fraction = dict_zones[n][z][2]
                                        price_material = price_material + fraction * pricelookup(prices, usage, material, zone, year, period, quantity ,true, true)
                                    end
                                    existings_price_per_kg_of_hydrogen[n, x, t] += price_material * quantity
                                end
                            end
                        end

                    end
                    end


                    ############################# EXPRESSIONS FOR CAPITAL COSTS ##################################
                    # Expressions for yearly capital costs based on new capacities added
                    println("CREATING EXPRESSIONS FOR CAPITAL AND FIXED COST")

                    @time begin
                    if i > 1
                        unregister(model, :yearly_cap_cost_production)
                        unregister(model, :yearly_fixed_cost_production)

                        unregister(model, :yearly_cap_cost_unextended)
                        unregister(model, :yearly_fixed_cost_unextended)

                        unregister(model, :yearly_cap_cost_storage)
                        unregister(model, :yearly_fixed_cost_storage)

                        unregister(model, :yearly_cap_cost_extended)
                        unregister(model, :yearly_fixed_cost_extended)

                    end


                    for n in node_set
                        for l in link_tot[n]
                            for p in pathway_set
                                for s in pathway_extended_storage_stages[p]
                                    for y in year_set

                                        for t in storage_constr_years
                                            if t >= y
                                                set_normalized_coefficient(total_storage_capacity_constr[t], pathway_extended_cap[n, l, p, s, y], normalized_ext_storage_cap[n, p, s, l, string(y)] * storage_available)
                                            end
                                        end

                                        for t in year_period_set
                                            #set_normalized_rhs(extended_storage_level_limit[n, l, p, s, y, t], available_extended_storage_cap[n, l, p, s, y, t[1]] * storage_available)
                                            set_normalized_coefficient(extended_storage_level_limit[n, l, p, s, y, t], pathway_extended_cap[n, l, p, s, y], - normalized_ext_storage_cap[n, p, s, l, string(y)] * storage_available)
                                        end
                                    end
                                end
                            end
                        end
                    end

                    JuMP.@expression(model, yearly_cap_cost_production[n in production_nodes, x in node_production_set[n], y in year_set, t in year_set], 0 + empty)
                    JuMP.@expression(model, yearly_fixed_cost_production[n in production_nodes, x in node_production_set[n], y in year_set, t in year_set], 0 + empty)

                    JuMP.@expression(model, yearly_cap_cost_unextended[n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set, t in year_set], 0 + empty)
                    JuMP.@expression(model, yearly_fixed_cost_unextended[n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set, t in year_set], 0 + empty)

                    JuMP.@expression(model, yearly_cap_cost_storage[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], y in year_set, t in year_set], 0 + empty)
                    JuMP.@expression(model, yearly_fixed_cost_storage[n in node_set, p in pathway_set, s in pathway_unextended_storage_stages[p], y in year_set, t in year_set], 0 + empty)

                    JuMP.@expression(model, yearly_cap_cost_extended[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set, t in year_set], 0 + empty)
                    JuMP.@expression(model, yearly_fixed_cost_extended[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set, t in year_set], 0 + empty)

                    for n in node_set
                        for y in year_set
                            if n in production_nodes
                                for x in node_production_set[n]
                                    life_time = production_data[i, n, x, string(y)].lifetime
                                    for t in y:min((y + life_time - 1), year_set[end])
                                        if t in year_set
                                            annualized_cap_cost = production_annualization[i, n, x, string(y)] * production_data[i, n, x, string(y)].capitalcost
                                            JuMP.add_to_expression!(yearly_cap_cost_production[n, x, y, t], new_production_cap[n, x, y] * annualized_cap_cost)
                                            JuMP.add_to_expression!(yearly_fixed_cost_production[n, x, y, t], new_production_cap[n, x, y] * production_data[i, n, x, string(y)].fixedcost)
                                        end
                                    end
                                end

                                for p in pathway_set
                                    for s in pathways_unextended[p]
                                        life_time = pathway_unext_flow_data[i, n, p, s, string(y)].lifetime
                                        for t in y:min((y + life_time - 1), year_set[end])
                                            if t in year_set
                                                annualized_cap_cost = pathway_unext_annualization[i, n, p, s, string(y)] * pathway_unext_flow_data[i, n, p, s, string(y)].capitalcost
                                                JuMP.add_to_expression!(yearly_cap_cost_unextended[n, p, s, y, t], pathway_unextended_cap[n, p, s, y] * annualized_cap_cost)
                                                JuMP.add_to_expression!(yearly_fixed_cost_unextended[n, p, s, y, t], pathway_unextended_cap[n, p, s, y] * pathway_unext_flow_data[i, n, p, s, string(y)].fixedcost)
                                            end
                                        end
                                    end


                                    for s in pathway_unextended_storage_stages[p]
                                        life_time = pathway_unext_storage_data[i, n, p, s, string(y)].lifetime
                                        for t in y:min((y + life_time - 1), year_set[end])
                                            if t in year_set
                                                annualized_cap_cost = pathway_storage_annualization[i, n, p, s, string(y)] * pathway_unext_storage_data[i, n, p, s, string(y)].capitalcost
                                                JuMP.add_to_expression!(yearly_cap_cost_storage[n, p, s, y, t], pathway_storage_cap[n, p, s, y] * annualized_cap_cost)
                                                JuMP.add_to_expression!(yearly_fixed_cost_storage[n, p, s, y, t], pathway_storage_cap[n, p, s, y] * pathway_unext_storage_data[i, n, p, s, string(y)].fixedcost)
                                            end
                                        end
                                    end

                                    for s in pathways_extended[p]
                                        for l in link_tot[n]
                                            life_time = pathway_ext_link_data[i, n, p, s, l, string(y)].lifetime
                                            for t in y:min((y + life_time - 1), year_set[end])
                                                if t in year_set
                                                    annualized_cap_cost = pathway_ext_annualization[i, n, p, s, l, string(y)] * pathway_ext_link_data[i, n, p, s, l, string(y)].capitalcost
                                                    JuMP.add_to_expression!(yearly_cap_cost_extended[n, l, p, s, y, t], pathway_extended_cap[n, l, p, s, y] * annualized_cap_cost)
                                                    JuMP.add_to_expression!(yearly_fixed_cost_extended[n, l, p, s, y, t], pathway_extended_cap[n, l, p, s, y] * pathway_ext_link_data[i, n, p, s, l, string(y)].fixedcost)
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    end
                    # Define objective function and solve optimization model
                    # Assumption - Storage IN and OUT costs are symmetric

                    println("BUILDING OBJECTIVE FUNCTION")
                    @time begin
                    o = create_objective_function(annualization,
                                    i,
                                    model,
                                    new_production_cap,pathway_unextended_cap, pathway_storage_cap, pathway_extended_cap,
                                    new_production, existing_production, pathway_stage_flow, pathway_stage_stor_in, pathway_stage_stor_out, link_flow,
                                    node_set, node_existings, production_nodes, link_tot, node_production_set, pathway_set,
                                    pathways_unextended, pathway_unextended_storage_stages,  pathways_extended,
                                    year_set, year_period_set, total_year_set, npv_array,
                                    yearly_cap_cost_production, yearly_cap_cost_unextended, yearly_cap_cost_storage, yearly_cap_cost_extended,
                                    yearly_fixed_cost_production, yearly_fixed_cost_unextended, yearly_fixed_cost_storage, yearly_fixed_cost_extended,
                                    production_data, production_price_per_kg_of_hydrogen,
                                    existings_data, existings_price_per_kg_of_hydrogen,
                                    pathway_unext_flow_data, pathway_unextended_price_per_kg_of_hydrogen,
                                    pathway_unext_storage_data, storage_price_per_kg_of_hydrogen,
                                    pathway_ext_link_data, pathway_extended_price_per_kg_of_hydrogen
                                                )
                    end

                    #println(o)
                    #println("Objective Function End")
                    println("Optimizing for iteration $(i)")
                    @time optimize!(model)
                    println(termination_status(model))
                    println(value(o))

                    map!(x -> 0.0, values(capacity_dict))

                    # Update capacities for each optimization iteration
                    for n in node_set
                        for y in year_set
                            for t in year_period_set
                                if n in production_nodes
                                    for x in node_production_set[n]
                                        infrastructureid = replace(string("INFR", '_', n, '_', x, '_', y), " " => "_")
                                        new_cap = round(value(new_production[n, x, y, t]))
                                        if new_cap > capacity_dict[infrastructureid]
                                            capacity_dict[infrastructureid] = new_cap
                                        end
                                    end
                                end

                                for p in pathway_set
                                    for s in pathways_unextended[p]
                                        infrastructureid = replace(string("INFR", '_', n, '_', p, '_', s, '_', y), " " => "_")
                                        new_cap = round(value(pathway_stage_flow[n, p, s, y, t]))
                                        if new_cap > capacity_dict[infrastructureid]
                                            capacity_dict[infrastructureid] = new_cap
                                        end
                                    end

                                    for s in pathways_extended[p]
                                        for l in link_tot[n]
                                            infrastructureid = replace(string("INFR", '_', n, '_', p,'_', s, '_', l, '_', y), " " => "_")
                                            new_cap = round(value(link_flow[n, l, p, s, y, t]))
                                            if new_cap > capacity_dict[infrastructureid]
                                                capacity_dict[infrastructureid] = new_cap
                                            end
                                        end
                                    end

                                    for s in pathway_unextended_storage_stages[p]
                                    infrastructureid = replace(string("INFR",  "_STORAGE_", n, '_', p, '_', s, '_', y), " " => "_")
                                    new_cap = round(value(pathway_stage_stor_in[n, p, s, y, t]))
                                            if new_cap > capacity_dict[infrastructureid]
                                                capacity_dict[infrastructureid] = new_cap
                                            end
                                    end
                                end
                            end
                        end

                    end

                    i += 1
                    current_objective = objective_value(model)
                    ϵ = abs((current_objective - prev_objective) / current_objective)
                    println("ϵ =  $(ϵ)")
                    prev_objective = current_objective
                end
                i -= 1
                println("ITERATIONS COMPLETED")

                map!(x->0, values(previous_ext_storage_cap))

                # Output Processing
                construction_output_year_set = filter(t -> t <= maximum(year_set), total_year_set)
                if last(year_set) == last_year
                    flow_output_year_set = year_set
                else
                    flow_output_year_set = sort(intersect(collect(current_year:min(total_year_set[end_index], current_year + rolling_step - 1)), total_year_set))
                end

                final_prices = DataFrame(Network_ID = String[], Year = Year[], Period = PeriodName[], Hydrogen_Price = Float64[])
                final_flow = DataFrame(Variable = String[], Infrastucture_ID = String[], Year = Year[], Period = PeriodName[], Technology = Technology[], Production_kg = Float64[], Flow_kg = Float64[], Loss_kg = Float64[], Cost = Float64[], Salvage_Value = Float64[])
                final_construction = DataFrame(Infrastucture_ID = String[], Network_ID = NetworkID[],  Technology = Technology[], Production = Productive[], Year = Year[], Lifetime = Int64[], Nameplate_Capacity = Float64[], Maximum_Utilization = Float64[], Length = Float64[], Capital_Cost = Float64[], Fixed_Operating_Cost = Float64[], Variable_Operating_Cost = Float64[], Storage_Technology = Union{String, Technology}[], Storage_Capacity = Float64[])

                println("WRITING RESULTS TO FILE")
                @time begin
                for n in demands_node_set
                    for t in flow_output_year_set
                        for h in period_set
                            if consumption[n, (t, h)] > 0.0
                                push!(final_prices, (n, t, h, round(value(dual(last_stage_to_consumption[n, (t, h)])), digits = 2)))
                            end
                        end
                    end
                end

                for t in flow_output_year_set
                    for n in node_set
                        for x in node_existings[n]
                            data = existings_data[n, x, (t, collect(keys(periods))[1])]
                            id = replace(string("INFR", "_", n, '_', x), " " => "_")
                            total_production = 0.0
                            total_cost = 0.0
                            for period in period_set
                                production = round(value(existing_production[n, x, (t, period)]))
                                cost = (data.cost + existings_price_per_kg_of_hydrogen[n, x, (t, period)]) * production
                                push!(final_flow, ("Existing", id, t, period, data.technology, production, 0.0, 0.0, cost, 0.0))
                                total_production += production
                                total_cost += cost
                            end

                            if t == start_year
                                push!(final_construction, (id, data.location, x, parseenum(Productive, "Central"), t[1], 1000, data.capacity, 1, 0, 0, 0, data.cost, "No Storage", 0.0))
                            end
                        end
                    end
                end

                for y in construction_output_year_set
                    for t in flow_output_year_set
                        for n in node_set
                            if n in production_nodes
                                for x in node_production_set[n]
                                        data = production_data[i, n, x, string(t)]
                                        id = replace(string("INFR", '_', n, '_', x, '_', y), " " => "_")
                                        cap_cost = 0.0

                                        yearly_fixed_cost = 0.0

                                        if y == t
                                            cap = value(new_production_cap[n, x, y])
                                            if cap > 1e-2
                                                cap_cost = data.capitalcost * cap
                                                fixed_cost_dict[id] = data.fixedcost * cap
                                                yearly_fixed_cost = fixed_cost_dict[id]
                                                push!(final_construction, (id, data.networkid, x, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, cap_cost, data.fixedcost * cap, data.variablecost, "No_Storage", 0.0))
                                                for m in t:min(t + data.lifetime - 1, total_year_set[end])
                                                    rolling_production_available_cap[n, x, string(t), string(m)] = cap
                                                end
                                            end

                                        end

                                        if in(id, keys(fixed_cost_dict)) && (t <= (y + data.lifetime - 1))
                                            total_cost = 0.0
                                            #total_cost = cap_cost + yearly_fixed_cost
                                            total_production = 0.0
                                            for period in period_set
                                                production = round(value(new_production[n, x, y, (t, period)]))
                                                cost = (data.variablecost + production_price_per_kg_of_hydrogen[n, x, (t, period)]) * production
                                                if production > 0.0
                                                    push!(final_flow, ("Production", id, t, period, data.technology, production, 0.0, 0.0, cost, 0.0))
                                                end
                                                total_production += production
                                                total_cost += cost
                                            end

                                        end

                                end
                            end

                            for p in pathway_set
                                for s in pathways_unextended[p]
                                    data = pathway_unext_flow_data[i, n, p, s, string(t)]
                                    id = replace(string("INFR", '_', n, '_', p, '_', s, '_', y), " " => "_")
                                    technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)

                                    cap_cost = 0.0
                                    yearly_fixed_cost = 0.0

                                    if y == t
                                        cap = value(pathway_unextended_cap[n, p, s, t])
                                        if cap > 1e-2
                                            storage = storagelookup(processlibrary, technology, t[1], 0.0, 0.0, id, n)
                                            cap_cost = data.capitalcost * cap
                                            fixed_cost_dict[id] = data.fixedcost * cap
                                            yearly_fixed_cost = fixed_cost_dict[id]
                                            push!(final_construction, (id, data.networkid, technology, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, data.capitalcost * cap, data.fixedcost * cap, data.variablecost, storage.technology, 0.0))

                                            for m in t:min(t + data.lifetime - 1, total_year_set[end])
                                                rolling_unext_available_cap[n, p, s, string(t), string(m)] = cap
                                            end
                                        end
                                    end


                                    if in(id, keys(fixed_cost_dict)) && (t <= (y + data.lifetime - 1))

                                        total_cost = 0.0

                                        #total_cost = cap_cost + yearly_fixed_cost
                                        total_flow = 0.0
                                        for period in period_set
                                            flow = round(value(pathway_stage_flow[n, p, s, y, (t, period)]))
                                            loss = flow * (1 - pathway_stage_yield[(p, s)])
                                            cost = (data.variablecost + pathway_unextended_price_per_kg_of_hydrogen[n, p, s, (t, period)]) * flow
                                            if flow > 0.0
                                                push!(final_flow, ("Unextended", id, t, period, data.technology, 0.0, flow, loss, cost, 0.0))
                                            end
                                            total_flow += flow
                                            total_cost += cost
                                        end

                                    end


                                end

                                for s in pathways_extended[p]
                                    for l in link_tot[n]
                                        data = pathway_ext_link_data[i, n, p, s, l, string(t)]
                                        id = replace(string("INFR", '_', n, '_', l, '_', p, '_', s, '_', y), " " => "_")
                                        technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)

                                        cap_cost = 0.0
                                        yearly_fixed_cost = 0.0

                                        if y == t
                                            cap = value(pathway_extended_cap[n, l, p, s, t])
                                            if cap > 1e-2
                                                storage = storagelookup(processlibrary, technology, t[1], 0.0, 0.0, id, n)
                                                cap_cost = data.capitalcost * cap
                                                fixed_cost_dict[id] = data.fixedcost * cap

                                                yearly_fixed_cost = fixed_cost_dict[id]
                                                if s in pathway_extended_storage_stages[p]
                                                    push!(final_construction, (id, data.networkid, technology, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, data.capitalcost * cap, data.fixedcost * cap, data.variablecost, storage.technology, normalized_ext_storage_cap[n, p, s, l, string(t)] * cap))
                                                else
                                                    push!(final_construction, (id, data.networkid, technology, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, data.capitalcost * cap, data.fixedcost * cap, data.variablecost, storage.technology, 0.0))
                                                end

                                                for m in t:min(t + data.lifetime - 1, total_year_set[end])
                                                    rolling_ext_available_cap[n, p, s, l, string(t), string(m)] = cap
                                                end

                                                if s in pathway_extended_storage_stages[p]
                                                    for m in t:min(t + data.lifetime - 1, total_year_set[end])
                                                        rolling_ext_storage_available_cap[n, p, s, l, string(t), string(m)] = normalized_ext_storage_cap[n, p, s, l, string(t)] * cap
                                                    end
                                                end
                                            end
                                        end


                                    if in(id, keys(fixed_cost_dict)) && (t <= (y + data.lifetime - 1))

                                        total_cost = 0.0
                                        #total_cost = cap_cost + yearly_fixed_cost
                                        total_flow = 0.0
                                        for period in period_set
                                            flow = round(value(link_flow[n, l, p, s, y, (t, period)]))
                                            cost = (data.variablecost + pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, (t, period)]) * flow
                                            loss = total_flow * (1 - pathway_stage_yield[(p, s)])
                                            if flow > 0.0
                                                push!(final_flow, ("Extended", id, t, period, data.technology, 0.0, flow, loss, cost, 0.0))
                                            end
                                            total_flow += flow
                                            total_cost += cost
                                        end


                                        if s in pathway_extended_storage_stages[p]
                                            for period in period_set
                                                stor_in = round(value(extended_stage_stor_in[n, l, p, s, y, (t, period)]))
                                                if stor_in > 0.0
                                                    push!(final_flow, ("Extended Storage In", id, t, period, data.technology, 0.0, stor_in, 0.0, 0.0, 0.0))
                                                end

                                                stor_out = round(value(extended_stage_stor_out[n, l, p, s, y, (t, period)]))
                                                if stor_out > 0.0
                                                    push!(final_flow, ("Extended Storage Out", id, t, period, data.technology, 0.0, stor_out, 0.0, 0.0, 0.0))
                                                end
                                            end

                                            if t == last(flow_output_year_set)
                                                rolling_extended_storage_level[n, l, p, s, string(y)] = round(value(extended_stage_stor_level[n, l, p, s, y, (t, final_period)]))
                                            end

                                        end
                                    end


                                    end
                                end

                                for s in pathway_unextended_storage_stages[p]
                                    data = pathway_unext_storage_data[i, n, p, s, string(t)]
                                    id = replace(string("INFR",  "_STORAGE_", n, '_', p, '_', s, '_', y), " " => "_")
                                    technology = getfield(processlibrary.pathways[PathwayKey(p, s)], :technology)

                                    cap_cost = 0.0

                                    if y == t
                                        cap = value(pathway_storage_cap[n, p, s, t])
                                        if cap > 1e-2
                                            storage = storagelookup(processlibrary, technology, t, 0.0, 0.0, id, n)
                                            cap_cost = data.capitalcost * cap
                                            fixed_cost_dict[id] = data.fixedcost * cap
                                            push!(final_construction, (id, data.networkid, storage.technology, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, data.capitalcost * cap, data.fixedcost * cap, data.variablecost, storage.technology, cap))

                                            for m in t:min(t + data.lifetime - 1, total_year_set[end])
                                                rolling_unext_storage_available_cap[n, p, s, string(t), string(m)] = cap
                                            end
                                        end
                                    end


                                    if in(id, keys(fixed_cost_dict)) && (t <= (y + data.lifetime - 1))

                                        #Storage In
                                        total_cost = 0.0
                                        #total_cost = cap_cost + yearly_fixed_cost
                                        total_flow = 0.0
                                        for period in period_set
                                            flow = round(value(pathway_stage_stor_in[n, p, s, y, (t, period)]))
                                            cost = (data.variablecost + storage_price_per_kg_of_hydrogen[n, p, s, (t, period)]) * flow
                                            loss = 0.0
                                            if flow > 0.0
                                                push!(final_flow, ("Storage In", id, t, period, data.technology, 0.0, flow, loss, cost, 0.0))
                                            end

                                            total_flow += flow
                                            total_cost += cost
                                        end

                                        #Storage Out
                                        total_cost = 0.0
                                        #total_cost = cap_cost + yearly_fixed_cost
                                        total_flow = 0.0
                                        for period in period_set
                                            flow = round(value(pathway_stage_stor_out[n, p, s, y, (t, period)]))
                                            cost = (data.variablecost + storage_price_per_kg_of_hydrogen[n, p, s, (t, period)]) * flow
                                            loss = 0.0
                                            if flow > 0.0
                                                push!(final_flow, ("Storage Out", id, t, period, data.technology, 0.0, flow, loss, cost, 0.0))
                                            end

                                            total_flow += flow
                                            total_cost += cost
                                        end

                                        if t == last(flow_output_year_set)
                                            rolling_storage_level[n, p, s, string(y)] = round(value(pathway_stage_stor_level[n, p, s, y, (t, final_period)]))
                                        end

                                    end

                                end

                            end

                        end
                    end
                end

                if isfile(outputs_flow)
                    stored_flow_df = CSV.read(outputs_flow, DataFrame)
                    flow_df_merge = vcat(stored_flow_df, final_flow)
                    CSV.write(outputs_flow, flow_df_merge, delim = "\t")

                    stored_construction_df = CSV.read(outputs_construction, DataFrame)
                    construction_df_merge = vcat(stored_construction_df, final_construction)
                    CSV.write(outputs_construction, construction_df_merge, delim = "\t")

                    stored_prices_df = CSV.read(outputs_prices, DataFrame)
                    prices_df_merge = vcat(stored_prices_df, final_prices)
                    CSV.write(outputs_prices, prices_df_merge, delim = "\t")


                else
                    CSV.write(outputs_flow, final_flow, delim = "\t")
                    CSV.write(outputs_construction, final_construction, delim = "\t")
                    CSV.write(outputs_prices, final_prices, delim = "\t")
                end
                end
            end
            model = nothing
            GC.gc()
        end
    end
    return
end



end
