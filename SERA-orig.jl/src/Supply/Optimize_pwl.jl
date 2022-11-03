module Optimize_pwl

using CSV
using DataFrames
using DataFramesMeta
using SERA.Supply.Types
using SERA.Supply.LookupFunctions
using SERA.Util
using JuMP
using GLPK
using Ipopt
using Xpress
using AxisArrays
using PiecewiseLinearOpt

const Xpress_optimizer = optimizer_with_attributes(Xpress.Optimizer,
                                                  "MIPRELSTOP" => 1e-3,
                                                  "BARGAPSTOP" => 1e-8,
                                                  "BARDUALSTOP" => 1e-8,
                                                  "BARPRIMALSTOP" => 1e-8,
                                                  "MATRIXTOL" => 1e-8,
                                                  "BARORDER" => 1,
                                                  "OUTPUTLOG" => 1)


import SERA.Supply.Types: MaterialName, MaterialUsage, Prices, ZoneID, Year, Network, ProcessLibrary, Demand, No, None, PeriodName, PathwayKey, Pathway, Construction

export makeModel_pwl

# Function that yields next period (season) in a given year or the first period of the next year

function next((year,period), period_set)
  if ((length(period_set) == 1) || (period == period_set[end]))
     result  = (year+1, period_set[1])
  else
     ind = findfirst(isequal(period), period_set)
     result = (year, period_set[ind+1])
  end
  return(result)
end

# Function that yields previous period (season) in a given year or the final period of the previous year

function prev((year,period), period_set)
  if ((length(period_set) == 1) || (period == period_set[1]))
     result  = (year-1, period_set[end])
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

function makeModel_pwl(network::Network,
    processlibrary::ProcessLibrary,
    demands::Demands, prices::Prices,
    usage::MaterialUsage, periods::Periods,
    niter::Int64, discount_rate::Float64,
    rollingWindow::Int64):: JuMP.Model

    model = Model(Xpress_optimizer)
    println("STARTING DATA PRE-PROCESSING")
    #model = Model(with_optimizer(Ipopt.Optimizer, max_iter=100))
    ################################## YEAR AND PERIOD SETS ##################################################

    # Create year, period and (year, period) sets
    demands_year_set = Vector{Year}()
    demands_period_set = Vector{PeriodName}()
    demands_year_period_set = Vector{Tuple{Year,PeriodName}}()

    for key in keys(demands)
	    push!(demands_year_set, key.year)
	    push!(demands_period_set, key.period)
        push!(demands_year_period_set, (key.year, key.period))
    end

    year_set = sort(unique(demands_year_set))
    period_set = sort(unique(demands_period_set))
    year_period_set = sort(unique(vcat(demands_year_period_set)))

    year_set_string = string.(year_set)
    start_year = year_set[1]

    npv_array = Dict(y => 1/(1 + discount_rate) ^ (y - start_year) for y in year_set)

    ################################## NODES ##################################################

    # Node numbers
    # All nodes in the network
    node_data = network.nodes
    node_set = collect(keys(node_data))

    # All production technologies in the network
    selected_costs = filter(kv->!(kv[2].productive in [No,None]), processlibrary.costs)
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
    demands_total = Dict{Tuple{NetworkID,Tuple{Year,PeriodName}},Float64}()

    for (key, value) in demands
        demands_total[(key.location, (key.year, key.period))] = (value.fuelConsumption + value.nonfuelConsumption) * periods[key.period]
    end

    ############################ PRODUCTION CAPABILITY ########################################
    production_nodes = filter(n -> node_data[n].productive != No && node_data[n].productive != None, node_set)

    ################################## EXISTINGS #############################################

    existings_capacity = Dict{Tuple{NetworkID, Technology, Tuple{Year,PeriodName}}, Float64}()
    existings_yield = Dict{Tuple{NetworkID, Technology, Tuple{Year,PeriodName}}, Float64}()
    existings_cost = Dict{Tuple{NetworkID, Technology, Tuple{Year,PeriodName}}, Float64}()

    for value in values(existings)
        existings_capacity[(value.location, value.technology, (value.year, value.period))] = value.capacity * periods[value.period]
        existings_yield[(value.location, value.technology, (value.year, value.period))] = value.yield
        existings_cost[(value.location, value.technology, (value.year, value.period))] = value.cost

        for year in year_set
            if (value.year !== year)
                for period in period_set
                    existings_capacity[(value.location, value.technology, (year, period))] = value.capacity
                    existings_yield[(value.location, value.technology, (year, period))] = value.yield
                    existings_cost[(value.location, value.technology, (year, period))] = value.cost
                end
            end
        end
    end

    ################################# LINKS ##################################################

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


    ################################# PATHWAYS ##################################################

    pathway_technology = Dict{Tuple{Pathway,Int64},Technology}()

    pathway_set = sort(unique(getfield.(keys(processlibrary.pathways), :pathway)))

    for (key,value) in processlibrary.pathways
	    pathway_technology[(key.pathway, key.stage)] = value.technology
    end

    pathway_stages = Dict{String, Int64}()
    pathway_storage_stages = Dict(p => Int64[] for p in pathway_set)
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
            storage = storagelookup(processlibrary, pathway_technology[pathway, stage], year_set[1], 1.0, 0.0, "", "")
            if (storage != "No Storage") && !isnothing(storage)
                push!(pathway_storage_stages[pathway], stage)
            end
        end

        for stage in stage_set
            pathway_stage_yield[(pathway, stage)] = getfield(processlibrary.pathways[PathwayKey(pathway, stage)], :yield)
            pathway_stage_yield[(pathway, stage)] = 1.0
        end
    end

    pathways_unextended_max_length = maximum(maximum.(collect(values(pathways_unextended))))
    pathways_extended_max_length = maximum(maximum.(collect(values(pathways_extended))))

    technology_storage = Dict{Technology, Union{Technology,Nothing}}()

    for (key,value) in processlibrary.costs
	      technology_storage[(key.technology)] = value.storage
    end

    production_data = AxisArray(Array{Construction}(undef, niter, length(node_set), length(production_set), length(year_set)), 1:niter, node_set, production_set, year_set_string)
    production_inputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(production_set), length(year_set)), 1:niter, node_set, production_set, year_set_string)
    production_outputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(production_set), length(year_set)), 1:niter, node_set, production_set, year_set_string)
    production_annualization = AxisArray(Array{Float64}(undef, niter, length(node_set), length(production_set), length(year_set)), 1:niter, node_set, production_set, year_set_string)

    existings_data = AxisArray(Array{Existing}(undef, niter, length(existings_node_set), length(existings_set), length(year_period_set)), 1:niter, existings_node_set, existings_set, year_period_set)
    existings_inputs = AxisArray(Array{Dict{}}(undef, niter, length(existings_node_set), length(existings_set), length(year_period_set)), 1:niter, existings_node_set, existings_set, year_period_set)
    existings_outputs = AxisArray(Array{Dict{}}(undef, niter, length(existings_node_set), length(existings_set), length(year_period_set)), 1:niter, existings_node_set, existings_set, year_period_set)

    pathway_unext_flow_data = AxisArray(Array{Construction}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)
    pathway_unext_flow_inputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)
    pathway_unext_flow_outputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)
    pathway_unext_annualization = AxisArray(Array{Float64}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)

    pathway_ext_link_data = AxisArray(Array{Construction}(undef, niter, length(node_set), length(pathway_set), pathways_extended_max_length, length(link_set), length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_extended_max_length, link_set, year_set_string)
    pathway_ext_link_inputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(pathway_set), pathways_extended_max_length, length(link_set), length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_extended_max_length, link_set, year_set_string)
    pathway_ext_link_outputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(pathway_set), pathways_extended_max_length, length(link_set), length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_extended_max_length, link_set, year_set_string)
    pathway_ext_annualization = AxisArray(Array{Float64}(undef, niter, length(node_set), length(pathway_set), pathways_extended_max_length, length(link_set), length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_extended_max_length, link_set, year_set_string)

    pathway_unext_storage_data = AxisArray(Array{Construction}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)
    pathway_storage_inputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)
    pathway_storage_outputs = AxisArray(Array{Dict{}}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)
    pathway_storage_annualization = AxisArray(Array{Float64}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)

    production_price_per_kg_of_hydrogen = AxisArray(Array{Float64}(undef, niter, length(node_set), length(production_set), length(year_set)), 1:niter, node_set, production_set, year_set_string)
    existings_price_per_kg_of_hydrogen = AxisArray(Array{Float64}(undef, niter, length(existings_node_set), length(existings_set), length(year_period_set)), 1:niter, existings_node_set, existings_set, year_period_set)
    pathway_unextended_price_per_kg_of_hydrogen =  AxisArray(Array{Float64}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)
    pathway_extended_price_per_kg_of_hydrogen = AxisArray(Array{Float64}(undef, niter, length(node_set), length(pathway_set), pathways_extended_max_length, length(link_set), length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_extended_max_length, link_set, year_set_string)
    storage_price_per_kg_of_hydrogen = AxisArray(Array{Float64}(undef, niter, length(node_set), length(pathway_set), pathways_unextended_max_length, length(year_set)), 1:niter, node_set, pathway_set, 1:pathways_unextended_max_length, year_set_string)

    ############################################# CREATE OPTIMIZATION MODEL ########################################################################


    println("STARTING OPTIMIZATION MODEL CREATION")
    ##### VARIABLES #################

    # Planning Variables:
    @variable(model, new_production_cap[n in production_nodes, x in production_set, t in year_set] >= 0.0)
    @variable(model, pathway_unextended_cap[n in node_set, p in pathway_set, s in pathways_unextended[p], t in year_set] >= 0.0)
    @variable(model, pathway_storage_cap[n in node_set, p in pathway_set, s in pathway_storage_stages[p], t in year_set] == 0.0)
    @variable(model, pathway_extended_cap[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], t in year_set] >= 0.0)

    # Operation Variables:
    # Production variable
    @variable(model, new_production[n in production_nodes, x in production_set, y in year_set, t in year_period_set] >= 0.0)

    # Existing variable
    @variable(model, existing_production[n in node_set, x in node_existings[n], t in year_period_set] >= 0.0)

    # Flow variables for all pathway stages
    @variable(model, pathway_stage_flow[n in node_set, p in pathway_set, s in 1:pathway_stages[p], y in year_set, t in year_period_set] >= 0.0)

    # Variables for hydrogen flowing into storage at unextended stages
    @variable(model, pathway_stage_stor_in[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_period_set] == 0.0)

    # Variables for hydrogen flowing out of storage at unextended stages
    @variable(model, pathway_stage_stor_out[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_period_set] >= 0.0)

    # Variables for hydrogen stored at unextended stages
    @variable(model, pathway_stage_stor_level[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_period_set] >= 0.0)

    # Variables for hydrogen flowing from node n across link l
    @variable(model, link_flow[n in node_set, link in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set, t in year_period_set] >= 0.0)

    # Dummy Variables
    @variable(model, empty == 0.0)

     ##### EXPRESSIONS #################

    # Expression that consumption equals total demands
    @expression(model, consumption[n in demands_node_set, t in year_period_set], demands_total[(n, t)])

    # Aggregate production expression
    @expression(model, aggregate_production[n in node_set, t in year_period_set], 1 * empty)

    # Expression for net hydrogen flow at each unextended stage of the pathway
    @expression(model, stage_net_flow[n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set, t in year_period_set], pathway_stage_yield[p, s] * pathway_stage_flow[n, p, s, y, t])

    # Expression for total available capacity based on new production added each year
    JuMP.@expression(model, available_production_cap[n in production_nodes, x in production_set, y in year_set, t in year_set], 0 + empty)

    # Expression for total available capacity based on new unextended pathways added each year
    JuMP.@expression(model, available_unextended_cap[n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set, t in year_set], 0 + empty)

    # Expression for total available capacity based on new storage pathways added each year
    JuMP.@expression(model, available_storage_cap[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_set], 0 + empty)

    # Expression for total available capacity based on new extended pathways added each year
    JuMP.@expression(model, available_extended_cap[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set, t in year_set], 0+ empty)


    println("POPULATING EXPRESSIONS")

    for n in node_set
        for t in year_period_set
            for x in node_existings[n]
                add_to_expression!(aggregate_production[n, t], existing_production[n, x, t])
            end
            if n in production_nodes
                add_to_expression!(aggregate_production[n, t], sum(new_production[n, x, y, t] for x in production_set, y in year_set))
            end
        end

        for p in pathway_set
            for s in pathway_storage_stages[p]
                for y in year_set
                    for t in year_period_set
                        add_to_expression!(stage_net_flow[n, p, s, y, t], pathway_stage_stor_out[n, p, s, y, t] - pathway_stage_stor_in[n, p, s, y, t])
                    end
                end
            end
         end

        for t in year_set
            for y in current_and_prev_years(t, year_set)
                if n in production_nodes
                    for x in production_set
                        if y + processcosts(processlibrary, x, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                            JuMP.add_to_expression!(available_production_cap[n, x, y, t], new_production_cap[n, x, y])
                        end
                    end
                end

                for p in pathway_set
                    for s in pathways_unextended[p]
                        dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                        technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                        if y + processcosts(processlibrary, technology, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                            JuMP.add_to_expression!(available_unextended_cap[n, p, s, y, t], pathway_unextended_cap[n, p, s, y])
                        end
                    end

                    for s in pathway_storage_stages[p]
                        dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                        technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                        storage = storagelookup(processlibrary, technology, y, 0.0, 0.0, "", n)
                        if y + processcosts(processlibrary, storage, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                            JuMP.add_to_expression!(available_storage_cap[n, p, s, y, t], pathway_storage_cap[n, p, s, y])
                        end
                    end

                    for s in pathways_extended[p]
                        for l in link_tot[n]
                            dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                            technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                            if y + processcosts(processlibrary, technology, y, 0.0, 0.0, "", n).lifetime - 1 >= t
                                JuMP.add_to_expression!(available_extended_cap[n, l, p, s, y, t], pathway_extended_cap[n, l, p, s, y])
                            end
                        end
                    end
                end
            end
        end
    end


    ##### CONSTRAINTS #################

    # Aggregate production constraint
    @constraint(model, aggregate_production_limit[n in node_set, t in year_period_set], aggregate_production[n, t] >= 0.0)

    # Sum of H2 flowing through the first stage of all pathways at node n should be equal to the total H2 production at node n
    @constraint(model, production_to_init_stage[n in node_set, t in year_period_set], aggregate_production[n, t] == sum(pathway_stage_flow[n, p, 1, y, t] for p in pathway_set, y in year_set))

    # H2 balance at unextended stages which are not the last stage of the pathway:
    @constraint(model, unext_flow_balance_with_stor[n in node_set, p in pathway_set, s in filter(x -> (x != last(pathways_unextended[p])), pathways_unextended[p]), y in year_set, t in year_period_set],
                stage_net_flow[n, p, s, y, t] == pathway_stage_flow[n, p, s + 1, y, t] )

    # H2 storage evolution at unextended stages (first time step):
    @constraint(model, unext_stor_evol_init[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set],
                pathway_stage_stor_level[n, p, s, y, year_period_set[1]] == pathway_stage_stor_in[n, p, s, y, year_period_set[1]] - pathway_stage_stor_out[n, p, s, y, year_period_set[1]] )

    # H2 storage evolution at unextended stages (subsequent time steps):
    @constraint(model, unext_stor_evol[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_period_set[2:end]],
                pathway_stage_stor_level[n, p, s, y, t] == pathway_stage_stor_level[n, p, s, y, prev(t, period_set)] + pathway_stage_stor_in[n, p, s, y, t] - pathway_stage_stor_out[n, p, s, y, t] )

    # H2 balance at extended stages:
    @constraint(model, ext_flow_balance[n in node_set, p in pathway_set, s in pathways_extended[p], y in year_set, t in year_period_set],
                pathway_stage_flow[n, p, s, y, t] + sum((pathway_stage_yield[p, s] * link_flow[other_node[(n, l)], l, p, s, y, t]) - link_flow[n, l, p, s, y, t] for l in link_tot[n])  == pathway_stage_flow[n, p, s + 1, y, t] )

    # Sum of net H2 flowing through the last stage of all pathways at node n should be equal to the total H2 demand at node n:
    @constraint(model, last_stage_to_consumption[n in demands_node_set, t in year_period_set],
                sum(stage_net_flow[n, p, pathway_stages[p], y, t] for p in pathway_set, y in year_set)  == consumption[n, t])


    # Capacity Bounds:
    # Constraint providing upper bound to production
	@constraint(model, production_constr[n in production_nodes, x in production_set, y in year_set, t in year_period_set], new_production[n, x, y, t] <= available_production_cap[n, x, y, t[1]] * periods[t[2]])

    # Constraint providing upper bound to existing
    @constraint(model, existings_constr[n in node_set, x in node_existings[n], t in year_period_set], existing_production[n, x, t] <= existings_capacity[(n, x, t)])

    # Constraint providing upper bound to pathway unextended stages
    @constraint(model, unext_cap_constr[n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set, t in year_period_set], pathway_stage_flow[n, p, s, y, t] <= available_unextended_cap[n, p, s, y, t[1]] * periods[t[2]])

    # Constraint providing upper bound to pathway storage stages
    @constraint(model, storage_cap_constr[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_period_set], pathway_stage_stor_level[n, p, s, y, t] <= available_storage_cap[n, p, s, y, t[1]])

    # Constraint providing upper bound to pathway extended stages
    @constraint(model, ext_cap_constr[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set, t in year_period_set], link_flow[n, l, p, s, y, t] <= available_extended_cap[n, l, p, s, y, t[1]] * periods[t[2]])

    println("POPULATED CORE OPTIMIZATION MODEL, STARTING ITERATIONS")


    # Objective Function of the model requires cost data which needs to be updated iteratively.

    # Initiate and populate data frame for flow variables before first iteration
    flowDF = DataFrame(Variable=String[], Infrastructure_ID=String[], Technology=Technology[], Pathway=Pathway[] , Stage=String[] , Node=NetworkID[], Link=NetworkID[], Year=Year[], Period=PeriodName[], Flow=Float64[], Operating_Cost=Float64[], Total_Cost=Float64[], Hydrogen_price=Float64[])

    for n in demands_node_set
        for t in year_period_set
            push!(flowDF, ("Consumption", "", "", "", "", n, "", t[1], t[2], 0.0, 0.0, 0.0, 0.0))
        end
    end

    infrastructure_id_set = String[]

    # Set initial flow for pathways
    init_flow = maximum(vcat(collect(values(demands_total)), collect(values(existings_capacity))))

    for n in node_set
        for y in year_set
            for t in year_period_set
                if n in production_nodes
                    for x in production_set
                        infrastructureid = replace(string("INFR", '_', n, '_', x, '_', y), " " => "_")
                        if !(infrastructureid in infrastructure_id_set)
                            push!(infrastructure_id_set, infrastructureid)
                        end
                        push!(flowDF, ("Production", infrastructureid, x, "", "", n, "", t[1], t[2], 0.0, 0.0, 0.0, 0.0))
                    end
                end

                for x in node_existings[n]
                    infrastructureid = replace(string("INFR", "_", n, '_', x), " " => "_")
                    if !(infrastructureid in infrastructure_id_set)
                        push!(infrastructure_id_set, infrastructureid)
                    end
                    push!(flowDF, ("Existing", infrastructureid, x, "", "", n, "", t[1], t[2], existings_capacity[n, x, t], existings_cost[n, x, t] * existings_capacity[n, x, t], 0.0, 0.0))
                end

                for p in pathway_set
                    for s in pathways_unextended[p]
                        infrastructureid = replace(string("INFR", '_', n, '_', p, '_', s, '_', y), " " => "_")
                        if !(infrastructureid in infrastructure_id_set)
                            push!(infrastructure_id_set, infrastructureid)
                        end
                        push!(flowDF, ("Unextended", infrastructureid, pathway_technology[(p, s)], p, string(s), n, "", t[1], t[2], init_flow, 0.0, 0.0, 0.0))
                    end

                    for s in pathways_extended[p]
                        for l in link_tot[n]
                            infrastructureid = replace(string("INFR", '_', n, '_', p,'_', s, '_', l, '_', y), " " => "_")
                            if !(infrastructureid in infrastructure_id_set)
                                push!(infrastructure_id_set, infrastructureid)
                            end
                            push!(flowDF, ("Link Flow", infrastructureid, pathway_technology[(p, s)], p, string(s), n, l, t[1], t[2], init_flow, 0.0, 0.0, 0.0))
                        end
                    end

                    for s in pathway_storage_stages[p]
                        infrastructureid = replace(string("INFR", "_STORAGE_", n, '_', p, '_', s, '_', y), " " => "_")
                        if !(infrastructureid in infrastructure_id_set)
                            push!(infrastructure_id_set, infrastructureid)
                        end
                        push!(flowDF, ("Storage In", infrastructureid, technology_storage[pathway_technology[(p, s)]], p, string(s), n, "", t[1], t[2], 0.0, 0.0, 0.0, 0.0))
                        push!(flowDF, ("Storage Out", infrastructureid, technology_storage[pathway_technology[(p, s)]], p, string(s), n, "", t[1], t[2], 0.0, 0.0, 0.0, 0.0))
                        push!(flowDF, ("Storage Level", infrastructureid, technology_storage[pathway_technology[(p, s)]], p, string(s), n, "", t[1], t[2], 0.0, 0.0, 0.0, 0.0))
                    end
                end
            end
        end
    end

    out_file = string("flowDF", 0,".csv")
    CSV.write(out_file, flowDF)
    niter = 1
    #Start Iterations
    for i in 1:niter
        println("READING COST DATA FOR ITERATION $(i)")
        for n in node_set
            for t in year_set
                if n in production_nodes
                    for x in production_set
                        infrastructureid = replace(string("INFR", '_', n, '_', x, '_', t), " " => "_")
                        sub1 = subdf_id(flowDF,infrastructureid)
                        #=
                        df_rows = Dict(y => filter(row -> (row[:Technology] == x && row[:Node] == n && row[:Year] ==y), flowDF) for y in year_set)
                        annual_flow = [sum(df_rows[y][:, :Flow]) for y in year_set]
                        capacity = maximum(annual_flow)
                        =#
                        capacity = maximum(flowDF[sub1, :Flow])
                        distance = 0.0
                        production_data[i, n, x, string(t)] = processcosts(processlibrary, x, t, capacity, distance, infrastructureid, n)

                        production_inputs[i, n, x, string(t)] = processinputs(processlibrary, x, t, capacity, distance)
                        production_outputs[i, n, x, string(t)] = processoutputs(processlibrary, x, t, capacity, distance)
                        input_materials = keys(production_inputs[i, n, x, string(t)])
                        production_price_per_kg_of_hydrogen[i, n, x, string(t)] = 0.0
                        for material in input_materials
                            # For this particular material, look up consumption/(kg of H2)
                            quantity = production_inputs[i, n, x, string(t)][material]
                            price_material = 0.0
                            for z in 1:length(dict_zones[n])
                                zone = dict_zones[n][z][1]
                                fraction = dict_zones[n][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            production_price_per_kg_of_hydrogen[i, n, x, string(t)] += price_material * quantity
                        end

                        output_materials = keys(production_outputs[i, n, x, string(t)])
                        for material in output_materials
                            # For this particular material, look up output/(kg of H2)
                            #println(material)
                            quantity = production_outputs[i, n, x, string(t)][material]
                            #println(quantity)
                            price_material = 0.0
                            for z in 1:length(dict_zones[n])
                                zone = dict_zones[n][z][1]
                                fraction = dict_zones[n][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            production_price_per_kg_of_hydrogen[i, n, x, string(t)] += price_material * quantity
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
                        sub1 = subdf_id(flowDF,infrastructureid)
                        sub2 = subdf2(flowDF, p, string(s))
                        capacity = maximum(flowDF[sub1 .& sub2, :Flow])/minimum(sort(unique(values(periods))))

                        distance = 0.0
                        dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                        technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                        pathway_unext_flow_data[i, n, p, s, string(t)] = processcosts(processlibrary, technology, t, capacity, distance, infrastructureid, n)

                        pathway_unext_flow_inputs[i, n, p, s, string(t)] = processinputs(processlibrary, technology, t, capacity, distance)
                        pathway_unext_flow_outputs[i, n, p, s, string(t)] = processoutputs(processlibrary, technology, t, capacity, distance)

                        input_materials = keys(pathway_unext_flow_inputs[i, n, p, s, string(t)])
                        pathway_unextended_price_per_kg_of_hydrogen[i, n, p, s, string(t)] = 0.0
                        for material in input_materials
                            # For this particular material, look up consumption/(kg of H2)
                            quantity = pathway_unext_flow_inputs[i, n, p, s, string(t)][material]
                            price_material = 0.0
                            for z in 1:length(dict_zones[n])
                                zone = dict_zones[n][z][1]
                                fraction = dict_zones[n][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            pathway_unextended_price_per_kg_of_hydrogen[i, n, p, s, string(t)] += price_material * quantity
                        end

                        output_materials = keys(pathway_unext_flow_outputs[i, n, p, s, string(t)])
                        for material in output_materials
                            # For this particular material, look up output/(kg of H2)
                            #println(material)
                            quantity = pathway_unext_flow_outputs[i, n, p, s, string(t)][material]
                            #println(quantity)
                            price_material = 0.0
                            for z in 1:length(dict_zones[n])
                                zone = dict_zones[n][z][1]
                                fraction = dict_zones[n][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            pathway_unextended_price_per_kg_of_hydrogen[i, n, p, s, string(t)] += price_material * quantity
                        end

                        if iszero(discount_rate)
                            pathway_unext_annualization[i, n, p, s, string(t)] = 1 / pathway_unext_flow_data[i, n, p, s, string(t)].lifetime
                        else
                            pathway_unext_annualization[i, n, p, s, string(t)] = discount_rate / (1 - (1 + discount_rate) ^ (-(pathway_unext_flow_data[i, n, p, s, string(t)].lifetime)))
                        end

                        storage = storagelookup(processlibrary, technology, t, capacity, distance, infrastructureid, n)
                        if ((storage != "No Storage") && !isnothing(storage))
                            infrastructureid = replace(string("INFR", "_STORAGE_", n, '_', p, '_', s, '_', t), " " => "_")
                            sub1 = subdf_id(flowDF,infrastructureid)
                            sub2 = subdf2(flowDF, p, string(s))
                          capacity = maximum(flowDF[sub1 .& sub2, :Flow])/minimum(sort(unique(values(periods))))
                          pathway_unext_storage_data[i, n, p, s, string(t)] = processcosts(processlibrary, storage, t, capacity, distance, infrastructureid, n)
                          pathway_storage_inputs[i, n, p, s, string(t)] = processinputs(processlibrary, storage, t, capacity, distance)
                          pathway_storage_outputs[i, n, p, s, string(t)] = processoutputs(processlibrary, storage, t, capacity, distance)

                          input_materials = keys(pathway_storage_inputs[i, n, p, s, string(t)])
                          storage_price_per_kg_of_hydrogen[i, n, p, s, string(t)] = 0.0
                          for material in input_materials
                          # For this particular material, look up consumption/(kg of H2)
                            quantity = pathway_storage_inputs[i, n, p, s, string(t)][material]
                            price_material = 0.0
                            for z in 1:length(dict_zones[n])
                                zone = dict_zones[n][z][1]
                                fraction = dict_zones[n][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            storage_price_per_kg_of_hydrogen[i, n, p, s, string(t)] += price_material * quantity
                          end

                          output_materials = keys(pathway_storage_outputs[i, n, p, s, string(t)])
                          for material in output_materials
                          # For this particular material, look up production/(kg of H2)
                            quantity = pathway_storage_outputs[i, n, p, s, string(t)][material]
                            price_material = 0.0
                            for z in 1:length(dict_zones[n])
                                zone = dict_zones[n][z][1]
                                fraction = dict_zones[n][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            storage_price_per_kg_of_hydrogen[i, n, p, s, string(t)] += price_material * quantity
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
                        capacity = maximum(filter(row -> (row[:Infrastructure_ID] == infrastructureid), flowDF)[:, :Flow])/minimum(sort(unique(values(periods))))
                        distance = link_length[l]

                        dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                        technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                        pathway_ext_link_data[i, n, p, s, l, string(t)] = processcosts(processlibrary, technology, t, capacity, distance, infrastructureid, l)

                        pathway_ext_link_inputs[i, n, p, s, l, string(t)] = processinputs(processlibrary, technology, t, capacity, distance)
                        pathway_ext_link_outputs[i, n, p, s, l, string(t)] = processoutputs(processlibrary, technology, t, capacity, distance)

                        input_materials = keys(pathway_ext_link_inputs[i, n, p, s, l, string(t)])
                        pathway_extended_price_per_kg_of_hydrogen[i, n, p, s, l, string(t)] = 0.0
                        for material in input_materials
                            # For this particular material, look up consumption/(kg of H2)
                            quantity = pathway_ext_link_inputs[i, n, p, s, l, string(t)][material]
                            price_material = 0.0
                            for z in 1:length(dict_zones[l])
                                zone = dict_zones[l][z][1]
                                fraction = dict_zones[l][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            pathway_extended_price_per_kg_of_hydrogen[i, n, p, s, l, string(t)] += price_material * quantity
                        end

                        output_materials = keys(pathway_ext_link_outputs[i, n, p, s, l, string(t)])
                        for material in output_materials
                            # For this particular material, look up production/(kg of H2)
                            quantity = pathway_ext_link_outputs[i, n, p, s, l, string(t)][material]
                            price_material = 0.0
                            for z in 1:length(dict_zones[l])
                                zone = dict_zones[l][z][1]
                                fraction = dict_zones[l][z][2]
                                price_material = price_material + fraction * pricelookup(prices, usage, material, zone, t, quantity ,true, true)
                            end
                            pathway_extended_price_per_kg_of_hydrogen[i, n, p, s, l, string(t)] += price_material * quantity
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

            for t in year_period_set
                for x in node_existings[n]
                    infrastructureid = replace(string("INFR", "_", n, '_', x), " " => "_")
                    sub1 = subdf_id(flowDF, infrastructureid)
                    capacity = existings_capacity[(n, x, t)] / periods[t[2]]
                    cost = existings_cost[n, x, t]
                    location = n
                    technology = x
                    year = t[1]
                    period = t[2]
                    distance = 0.0
                    yield = existings_yield[n, x, t]
                    existings_data[i, n, x, t] = Existing(location, technology, year, period, capacity, yield, cost)

                    existings_inputs[i, n, x, t] = processinputs(processlibrary, x, year, capacity, distance)
                    existings_outputs[i, n, x, t] = processoutputs(processlibrary, x, year, capacity, distance)
                    input_materials = keys(existings_inputs[i, n, x, t])
                    existings_price_per_kg_of_hydrogen[i, n, x, t] = 0.0
                    for material in input_materials
                        # For this particular material, look up consumption/(kg of H2)
                          quantity = existings_inputs[i, n, x, t][material]
                          price_material = 0.0
                          for z in 1:length(dict_zones[n])
                              zone = dict_zones[n][z][1]
                              fraction = dict_zones[n][z][2]
                              price_material = price_material + fraction * pricelookup(prices, usage, material, zone, year, quantity ,true, true)
                          end
                          existings_price_per_kg_of_hydrogen[i, n, x, t] += price_material * quantity
                    end

                    output_materials = keys(existings_outputs[i, n, x, t])
                    for material in output_materials
                        # For this particular material, look up production/(kg of H2)
                          quantity = existings_outputs[i, n, x, t][material]
                          price_material = 0.0
                          for z in 1:length(dict_zones[n])
                              zone = dict_zones[n][z][1]
                              fraction = dict_zones[n][z][2]
                              price_material = price_material + fraction * pricelookup(prices, usage, material, zone, year, quantity ,true, true)
                          end
                          existings_price_per_kg_of_hydrogen[i, n, x, t] += price_material * quantity
                    end
                end
            end
        end

        ############################# EXPRESSIONS FOR CAPITAL COSTS ##################################
        # Expressions for yearly capital costs based on new capacities added

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

        JuMP.@expression(model, yearly_cap_cost_production[n in production_nodes, x in production_set, y in year_set, t in year_set], 0 + empty)
        JuMP.@expression(model, yearly_fixed_cost_production[n in production_nodes, x in production_set, y in year_set, t in year_set], 0 + empty)

        JuMP.@expression(model, yearly_cap_cost_unextended[n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set, t in year_set], 0 + empty)
        JuMP.@expression(model, yearly_fixed_cost_unextended[n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set, t in year_set], 0 + empty)

        JuMP.@expression(model, yearly_cap_cost_storage[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_set], 0 + empty)
        JuMP.@expression(model, yearly_fixed_cost_storage[n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set, t in year_set], 0 + empty)

        JuMP.@expression(model, yearly_cap_cost_extended[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set, t in year_set], 0 + empty)
        JuMP.@expression(model, yearly_fixed_cost_extended[n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set, t in year_set], 0 + empty)


        for n in node_set
            for y in year_set
                if n in production_nodes
                    for x in production_set
                        life_time = production_data[i, n, x, string(y)].lifetime
                        for t in y:min((y + life_time - 1), year_set[end])
                            annualized_cap_cost = production_annualization[i, n, x, string(y)] * production_data[i, n, x, string(y)].capitalcost
                            JuMP.add_to_expression!(yearly_cap_cost_production[n, x, y, t], new_production_cap[n, x, y] * annualized_cap_cost)
                            JuMP.add_to_expression!(yearly_fixed_cost_production[n, x, y, t], new_production_cap[n, x, y] * production_data[i, n, x, string(y)].fixedcost)
                        end
                    end

                    for p in pathway_set
                        for s in pathways_unextended[p]
                            life_time = pathway_unext_flow_data[i, n, p, s, string(y)].lifetime
                            for t in y:min((y + life_time - 1), year_set[end])
                                annualized_cap_cost = pathway_unext_annualization[i, n, p, s, string(y)] * pathway_unext_flow_data[i, n, p, s, string(y)].capitalcost
                                JuMP.add_to_expression!(yearly_cap_cost_unextended[n, p, s, y, t], pathway_unextended_cap[n, p, s, y] * annualized_cap_cost)
                                JuMP.add_to_expression!(yearly_fixed_cost_unextended[n, p, s, y, t], pathway_unextended_cap[n, p, s, y] * pathway_unext_flow_data[i, n, p, s, string(y)].fixedcost)
                            end
                        end

                        for s in pathway_storage_stages[p]
                            life_time = pathway_unext_storage_data[i, n, p, s, string(y)].lifetime
                            for t in y:min((y + life_time - 1), year_set[end])
                                annualized_cap_cost = pathway_storage_annualization[i, n, p, s, string(y)] * pathway_unext_storage_data[i, n, p, s, string(y)].capitalcost
                                JuMP.add_to_expression!(yearly_cap_cost_storage[n, p, s, y, t], pathway_storage_cap[n, p, s, y] * annualized_cap_cost)
                                JuMP.add_to_expression!(yearly_fixed_cost_storage[n, p, s, y, t], pathway_storage_cap[n, p, s, y] * pathway_unext_storage_data[i, n, p, s, string(y)].fixedcost)
                            end
                        end

                        for s in pathways_extended[p]
                            for l in link_tot[n]
                                life_time = pathway_ext_link_data[i, n, p, s, l, string(y)].lifetime
                                for t in y:min((y + life_time - 1), year_set[end])
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

        # Define objective function and solve optimization model
        # Assumption - Storage IN and OUT costs are symmetric

        println("Solve model")
        o=@objective(model, Min,
        #Capital Costs
                                  sum(npv_array[t] * (
                                    sum(piecewiselinear(model, new_production_cap[n, x, t], production_data[i, n, x, string(t)].pwl_cap_cost[1], production_data[i, n, x, string(t)].pwl_cap_cost[2]) for n in production_nodes, x in production_set) +
                                    sum(piecewiselinear(model, pathway_unextended_cap[n, p, s, t], pathway_unext_flow_data[i, n, p, s, string(t)].pwl_cap_cost[1], pathway_unext_flow_data[i, n, p, s, string(t)].pwl_cap_cost[2]) for n in node_set, p in pathway_set, s in pathways_unextended[p]) +
                                    sum(piecewiselinear(model, pathway_storage_cap[n, p, s, t], pathway_unext_storage_data[i, n, p, s, string(t)].pwl_cap_cost[1], pathway_unext_storage_data[i, n, p, s, string(t)].pwl_cap_cost[2]) for n in node_set, p in pathway_set, s in pathway_storage_stages[p]) +
                                    sum(piecewiselinear(model, pathway_extended_cap[n, l, p, s, t], pathway_ext_link_data[i, n, p, s, l, string(t)].pwl_cap_cost[1], pathway_ext_link_data[i, n, p, s, l, string(t)].pwl_cap_cost[2]) for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p])
                                  )
                                    for t in year_set)+

        #Fixed Costs
                                    sum(npv_array[t] * (
                                        sum(piecewiselinear(model, available_production_cap[n, x, y, t], production_data[i, n, x, string(t)].pwl_fix_cost[1], production_data[i, n, x, string(t)].pwl_fix_cost[2]) for n in production_nodes, x in production_set) +
                                        sum(piecewiselinear(model, available_unextended_cap[n, p, s, y, t], pathway_unext_flow_data[i, n, p, s, string(t)].pwl_fix_cost[1], pathway_unext_flow_data[i, n, p, s, string(t)].pwl_fix_cost[2]) for n in node_set, p in pathway_set, s in pathways_unextended[p]) +
                                        sum(piecewiselinear(model, available_storage_cap[n, p, s, y, t], pathway_unext_storage_data[i, n, p, s, string(t)].pwl_fix_cost[1], pathway_unext_storage_data[i, n, p, s, string(t)].pwl_fix_cost[2]) for n in node_set, p in pathway_set, s in pathway_storage_stages[p]) +
                                        sum(piecewiselinear(model, available_extended_cap[n, l, p, s, y, t], pathway_ext_link_data[i, n, p, s, l, string(t)].pwl_fix_cost[1], pathway_ext_link_data[i, n, p, s, l, string(t)].pwl_fix_cost[2]) for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p])
                                                    )
                                  for y in year_set, t in year_set) +
        #Operation Costs
                                  sum(npv_array[t[1]] * (
                                    sum((production_data[i, n, x, string(t[1])].variablecost + production_price_per_kg_of_hydrogen[i, n, x, string(t[1])]) * new_production[n, x, y, t] for n in production_nodes, x in production_set, y in year_set) +
                                    sum((existings_data[i, n, x, t].cost + existings_price_per_kg_of_hydrogen[i, n, x, t]) * existing_production[n, x, t] for n in node_set, x in node_existings[n]) +
                                    sum((pathway_unext_flow_data[i, n, p, s, string(t[1])].variablecost + pathway_unextended_price_per_kg_of_hydrogen[i, n, p, s, string(t[1])]) * pathway_stage_flow[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathways_unextended[p], y in year_set) +
                                    sum((pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[i, n, p, s, string(t[1])]) * pathway_stage_stor_in[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set) +
                                    sum((pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[i, n, p, s, string(t[1])]) * pathway_stage_stor_out[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in year_set) +
                                    sum((pathway_ext_link_data[i, n, p, s, l, string(t[1])].variablecost + pathway_extended_price_per_kg_of_hydrogen[i, n, p, s, l, string(t[1])]) * link_flow[n, l, p, s, y, t]
                                        for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in year_set)
                                                        )
                                  for t in year_period_set))
        #println("Objective Function Start")
        #println(o)
        #println("Objective Function End")
        println("Optimizing for iteration $(i)")
        @time optimize!(model)
        println(termination_status(model))


        # Update data frame for each optimization iteration

        for n in demands_node_set
          for t in year_period_set
            sub0 = subdf0(flowDF,"Consumption", n ,t[1], t[2])
            flowDF[sub0, :Flow] .= round(value(consumption[n, t]))
            #flowDF[sub0, :Hydrogen_price] .= dual(last_stage_to_consumption[n, t])
          end
        end


        for n in node_set
            for t in year_period_set
                for x in node_existings[n]
                    sub1 = subdf1(flowDF, "Existing", x, n, t[1])
                    sub2 = (flowDF.Period .== t[2])
                    flowDF[sub1 .& sub2, :Flow] .= round(value(existing_production[n, x, t]))
                    flowDF[sub1 .& sub2, :Operating_Cost] .= (existings_data[i, n, x, t].cost + existings_price_per_kg_of_hydrogen[i, n, x, t]) * flowDF[sub1 .& sub2, :Flow]
                    flowDF[sub1 .& sub2, :Total_Cost] .= flowDF[sub1 .& sub2, :Operating_Cost]
                end
            end

            for y in year_set
                for t in year_period_set
                    if n in production_nodes
                        for x in production_set
                            infrastructureid = replace(string("INFR", '_', n, '_', x, '_', y), " " => "_")
                            sub1 = subdf_id(flowDF, infrastructureid)
                            sub2 = subdf4(flowDF, t)
                            flowDF[sub1 .& sub2, :Flow] .= round(value(new_production[n, x, y, t]))
                            flowDF[sub1 .& sub2, :Operating_Cost] .= (production_data[i, n, x, string(t[1])].variablecost + production_price_per_kg_of_hydrogen[i, n, x, string(t[1])]) * flowDF[sub1 .& sub2, :Flow]

                            if t[1] == y
                                #cap_cost_add = production_data[i, n, x, string(y)].capitalcost * value(new_production_cap[n, x, y])
                                cap_cost_add = extrap_cost(value(new_production_cap[n, x, y]), production_data[i, n, x, string(y)].pwl_cap_cost[1], production_data[i, n, x, string(y)].pwl_cap_cost[2])
                            else
                                cap_cost_add = 0.0
                            end

                            flowDF[sub1 .& sub2, :Total_Cost] .= flowDF[sub1 .& sub2, :Operating_Cost] + [cap_cost_add + value(yearly_fixed_cost_production[n, x, y, t[1]])]
                        end
                    end

                    for p in pathway_set
                        for s in pathways_unextended[p]
                            infrastructureid = replace(string("INFR", '_', n, '_', p, '_', s, '_', y), " " => "_")
                            sub1 = subdf_id(flowDF, infrastructureid)
                            sub2 = subdf4(flowDF, t)
                            flowDF[sub1 .& sub2, :Flow] .= round(value(pathway_stage_flow[n, p, s, y, t]))

                            if t[1] == y
                                #cap_cost_add = pathway_unext_flow_data[i, n, p, s, string(y)].capitalcost * value(pathway_unextended_cap[n, p, s, y])
                                cap_cost_add = extrap_cost(value(pathway_unextended_cap[n, p, s, y]), pathway_unext_flow_data[i, n, p, s, string(y)].pwl_cap_cost[1], pathway_unext_flow_data[i, n, p, s, string(y)].pwl_cap_cost[2])
                            else
                                cap_cost_add = 0.0
                            end

                            flowDF[sub1 .& sub2, :Operating_Cost] .= (pathway_unext_flow_data[i, n, p, s, string(t[1])].variablecost + pathway_unextended_price_per_kg_of_hydrogen[i, n, p, s, string(t[1])]) * flowDF[sub1 .& sub2, :Flow]
                            flowDF[sub1 .& sub2, :Total_Cost] .= flowDF[sub1 .& sub2, :Operating_Cost] + [cap_cost_add + value(yearly_fixed_cost_unextended[n, p, s, y, t[1]])]
                        end

                        for s in pathways_extended[p]
                            for l in link_tot[n]
                            infrastructureid = replace(string("INFR", '_', n, '_', p,'_', s, '_', l, '_', y), " " => "_")
                            sub1 = subdf_id(flowDF, infrastructureid)
                            sub2 = subdf4(flowDF, t)
                            flowDF[sub1 .& sub2, :Flow] .= round(value(link_flow[n, l, p, s, y, t]))

                            if t[1] == y
                                #cap_cost_add = pathway_ext_link_data[i, n, p, s, l, string(y)].capitalcost * value(pathway_extended_cap[n, l, p, s, y])
                                cap_cost_add = extrap_cost(value(pathway_extended_cap[n, l, p, s, y]), pathway_ext_link_data[i, n, p, s, l, string(y)].pwl_cap_cost[1], pathway_ext_link_data[i, n, p, s, l, string(y)].pwl_cap_cost[2])
                            else
                                cap_cost_add = 0.0
                            end

                            flowDF[sub1 .& sub2, :Operating_Cost] .= (pathway_ext_link_data[i, n, p, s, l, string(t[1])].variablecost+ pathway_extended_price_per_kg_of_hydrogen[i, n, p, s, l, string(t[1])]) * flowDF[sub1 .& sub2, :Flow]
                            flowDF[sub1 .& sub2, :Total_Cost] .= flowDF[sub1 .& sub2, :Operating_Cost] + [cap_cost_add + value(yearly_fixed_cost_extended[n, l, p, s, y, t[1]])]
                            end
                        end

                        for s in pathway_storage_stages[p]
                        infrastructureid = replace(string("INFR",  "_STORAGE_", n, '_', p, '_', s, '_', y), " " => "_")
                        sub1 = subdf_id(flowDF, infrastructureid)
                        sub2 = subdf4(flowDF, t)
                        sub3 = subdf5(flowDF, "Storage In")
                        flowDF[sub1 .& sub2 .& sub3, :Flow] .= round(value(pathway_stage_stor_in[n, p, s, y, t]))

                        if t[1] == y
                            #cap_cost_add = pathway_unext_storage_data[i, n, p, s, string(y)].capitalcost * value(pathway_storage_cap[n, p, s, y])
                            cap_cost_add = extrap_cost(value(pathway_storage_cap[n, p, s, y]), pathway_unext_storage_data[i, n, p, s, string(y)].pwl_cap_cost[1], pathway_unext_storage_data[i, n, p, s, string(y)].pwl_cap_cost[2])
                        else
                            cap_cost_add = 0.0
                        end

                        flowDF[sub1 .& sub2 .& sub3, :Operating_Cost] .= (pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[i, n, p, s, string(t[1])]) * flowDF[sub1 .& sub2 .& sub3, :Flow]
                        flowDF[sub1 .& sub2 .& sub3, :Total_Cost] .= flowDF[sub1 .& sub2 .& sub3, :Operating_Cost] + [cap_cost_add + value(yearly_fixed_cost_storage[n, p, s, y, t[1]])]

                        sub3 = subdf5(flowDF, "Storage Out")
                        flowDF[sub1 .& sub2 .& sub3, :Flow] .= round(value(pathway_stage_stor_out[n, p, s, y, t]))
                        flowDF[sub1 .& sub2 .& sub3, :Operating_Cost] .= (pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[i, n, p, s, string(t[1])]) * flowDF[sub1 .& sub2 .& sub3, :Flow]

                        sub3 = subdf5(flowDF, "Storage Level")
                        flowDF[sub1 .& sub2 .& sub3, :Flow] .= round(value(pathway_stage_stor_level[n, p, s, y, t]))
                        end
                    end
                end
            end
        end

        # Output data frame for each iteration

        out_file = string("flowDF",i,".csv")
        CSV.write(out_file, flowDF)
    end

    # Output Processing
    final_flow = DataFrame(Variable = String[], Infrastucture_ID = String[], Year = Year[], Technology = Technology[], Production_kg = Float64[], Flow_kg = Float64[], Loss_kg = Float64[], Cost = Float64[], Salvage_Value = Float64[])

    final_construction = DataFrame(Infrastucture_ID = String[], Network_ID = NetworkID[],  Technology = Technology[], Production = Productive[], Year = Year[], Lifetime = Int64[], Nameplate_Capacity = Float64[], Maximum_Utilization = Float64[], Length = Float64[], Capital_Cost = Float64[], Fixed_Operating_Cost = Float64[], Variable_Operating_Cost = Float64[], Storage_Technology = Union{String, Technology}[])

    variable_set = ["Production", "Existing", "Unextended", "Link Flow", "Storage In", "Storage Out", "Storage Level"]

    for id in infrastructure_id_set
        for t in year_set
            for v in variable_set
                data = filter(row -> (row[:Infrastructure_ID] == id && row[:Year] == t && row[:Variable] == v), flowDF)
                unique!(data)
                if nrow(data) > 0
                    total_flow = 0.0
                    total_cost = 0.0
                    for r in 1:nrow(data)
                        total_flow += data[r, :Flow]
                        if r == 1
                            total_cost += data[r, :Total_Cost]
                        else
                            total_cost += data[r, :Operating_Cost]
                        end
                    end
                    if total_flow > 0 || data[1, :Variable] == "Existing"
                        if data[1, :Variable] == "Production" || data[1, :Variable] == "Existing"
                            production = total_flow
                            flow = 0.0
                            loss = 0.0
                        elseif data[1, :Variable] == "Unextended" || data[1, :Variable] == "Link Flow" || data[1, :Variable] == "Storage In" || data[1, :Variable] == "Storage Out" || data[1, :Variable] == "Storage Level"
                            production = 0.0
                            flow = total_flow
                            loss = (1 - pathway_stage_yield[(data[1, :Pathway], parse(Int64, data[1, :Stage]))]) * flow
                        end
                        push!(final_flow, (v, id, t, data[1, :Technology], production, flow, loss, total_cost, 0.0))
                    end
                end
            end
        end
    end

    for t in year_set
        for n in node_set
            for x in node_existings[n]
                if t == start_year
                    data = existings_data[niter, n, x, (t, collect(keys(periods))[1])]
                    id = replace(string("INFR", "_", n, '_', x), " " => "_")
                    push!(final_construction, (id, data.location, x, parseenum(Productive, "Central"), t[1], 1000, data.capacity, 1, 0, 0, 0, data.cost, "No Storage"))
                end
            end

            if n in production_nodes
                for x in production_set
                    cap = value(new_production_cap[n, x, t])
                    if cap > 1e-4
                        data = production_data[niter, n, x, string(t)]
                        id = replace(string("INFR", '_', n, '_', x, '_', t), " " => "_")
                        cap_cost = extrap_cost(cap, data.pwl_cap_cost[1], data.pwl_cap_cost[2])
                        fix_cost = extrap_cost(cap, data.pwl_fix_cost[1], data.pwl_fix_cost[2])
                        #push!(final_construction, (id, data.networkid, x, data.productive, t, data.lifetime, value(new_production_cap[n, x, t]), data.dutycycle, data.length, data.capitalcost * value(new_production_cap[n, x, t]), data.fixedcost * value(new_production_cap[n, x, t]), data.variablecost, "No_Storage"))
                        push!(final_construction, (id, data.networkid, x, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, cap_cost, fix_cost, data.variablecost, "No_Storage"))
                    end
                end
            end

            for p in pathway_set
                for s in pathways_unextended[p]
                    cap = value(pathway_unextended_cap[n, p, s, t])
                    if cap > 1e-4
                        data = pathway_unext_flow_data[niter, n, p, s, string(t)]
                        id = replace(string("INFR", '_', n, '_', p, '_', s, '_', t), " " => "_")
                        dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                        technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                        storage = storagelookup(processlibrary, technology, t[1], 0.0, 0.0, id, n)
                        cap_cost = extrap_cost(cap, data.pwl_cap_cost[1], data.pwl_cap_cost[2])
                        fix_cost = extrap_cost(cap, data.pwl_fix_cost[1], data.pwl_fix_cost[2])
                        #push!(final_construction, (id, data.networkid, technology, data.productive, t, data.lifetime, value(pathway_unextended_cap[n, p, s, t]), data.dutycycle, data.length, data.capitalcost * value(pathway_unextended_cap[n, p, s, t]), data.fixedcost * value(pathway_unextended_cap[n, p, s, t]), data.variablecost, storage))
                        push!(final_construction, (id, data.networkid, technology, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, cap_cost, fix_cost, data.variablecost, storage))
                    end
                end

                for s in pathways_extended[p]
                    for l in link_tot[n]
                        cap = value(pathway_extended_cap[n, l, p, s, t])
                        if cap > 1e-4
                            data = pathway_ext_link_data[niter, n, p, s, l, string(t)]
                            id = replace(string("INFR", '_', n, '_', l, '_', p, '_', s, '_', t), " " => "_")
                            dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                            technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                            storage = storagelookup(processlibrary, technology, t, 0.0, 0.0, id, n)
                            cap_cost = extrap_cost(cap, data.pwl_cap_cost[1], data.pwl_cap_cost[2])
                            fix_cost = extrap_cost(cap, data.pwl_fix_cost[1], data.pwl_fix_cost[2])
                            #push!(final_construction, (id, data.networkid, technology, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, data.capitalcost * cap, data.fixedcost * cap, data.variablecost, storage))
                            push!(final_construction, (id, data.networkid, technology, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, cap_cost, fix_cost, data.variablecost, storage))
                        end
                    end
                end

                for s in pathway_storage_stages[p]
                    cap = value(pathway_storage_cap[n, p, s, t])
                    if cap > 1e-4
                        data = pathway_unext_storage_data[niter, n, p, s, string(t)]
                        id = replace(string("INFR",  "_STORAGE_", n, '_', p, '_', s, '_', t), " " => "_")
                        dict_path_spec = filter(kv ->(kv[1].pathway == p && kv[1].stage == s), processlibrary.pathways)
                        technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
                        storage = storagelookup(processlibrary, technology, t, 0.0, 0.0, id, n)
                        cap_cost = extrap_cost(cap, data.pwl_cap_cost[1], data.pwl_cap_cost[2])
                        fix_cost = extrap_cost(cap, data.pwl_fix_cost[1], data.pwl_fix_cost[2])
                        #push!(final_construction, (id, data.networkid, storage, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, data.capitalcost * cap, data.fixedcost * cap, data.variablecost, storage))
                        push!(final_construction, (id, data.networkid, storage, data.productive, t, data.lifetime, cap, data.dutycycle, data.length, cap_cost, fix_cost, data.variablecost, storage))
                    end
                end

            end
        end
    end

    outputs_flow = joinpath("outputs", "flow_SERA_2_pwl.tsv")
    CSV.write(outputs_flow, final_flow, delim = "\t")



    outputs_construction = joinpath("outputs", "construction_SERA_2_pwl.tsv")
    CSV.write(outputs_construction, final_construction, delim = "\t")

    return model
end



end
