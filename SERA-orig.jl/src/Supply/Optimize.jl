module Optimize

using CSV
using DataFrames
using DataFramesMeta
using SERA.Supply.Types
using SERA.Supply.LookupFunctions
using SERA.Util
using JuMP
using GLPK
using Ipopt

import SERA.Supply.Types: MaterialName, MaterialUsage, Prices, ZoneID, Year, Network, ProcessLibrary, Demand, No, None, PeriodName, PathwayKey, Pathway, Construction

export makeModel

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

# Formulate Optimization Model

function makeModel(network::Network, processlibrary::ProcessLibrary, demands::Demands, prices::Prices, usage::MaterialUsage, periods::Periods, niter::Int64):: JuMP.Model
    model = Model(Ipopt.Optimizer)
    #model = Model(with_optimizer(Ipopt.Optimizer, max_iter=100))

    # Node numbers
    # All nodes in the network
    node_set = collect(keys(network.nodes))
    selected_costs = filter(kv->!(kv[2].productive in [No,None]), processlibrary.costs)

    # All production technologies in the network
    production_set = sort(unique(getfield.(keys(selected_costs), :technology)))

    # All existing technologies in the network
    existings = network.existings

    println("Existings")
    println(existings)

    existings_set = sort(unique(getfield.(values(existings), :technology)))

    # All nodes with existing capacity
    existings_node_set = sort(unique(getfield.(values(existings), :location)))
    # All nodes with demands
    demands_node_set = sort(unique(getfield.(keys(demands), :location)))

    println("Production Set")
    println(production_set)

    println("Node Set")
    println(node_set)

    println("Existings Set")
    println(existings_set)

    println("Existings Node Set")
    println(existings_node_set)

    println("Demands Node Set")
    println(demands_node_set)

    # Dictionary for links going to a node
    link_in = Dict{}()
    # Dictionary for links leaving a node
    link_out = Dict{}()
    # Dictionary for all links associated with a node
    link_tot = Dict{}()
    link_length = Dict{}()
    # Dictionary to find node at other end of link from given node
    other_node = Dict{}()
    # Dictionary connecting nodes to zones
    dict_zones = Dict{}()

    # Create set of network IDs for links
    link_set = collect(keys(network.links))

    println("Link Lengths")
    for link in link_set
        link_length[link] = network.links[link].length
        println(link,",",link_length[link])
    end

    println("Network Links")
    println(network.links)
    # Populate the link related dictionaries
    for node in node_set
        dict_links_into_node = filter(kv -> (kv[2].to == node), network.links)
        # Links going to a node
        links_into_node = sort(unique(keys(dict_links_into_node)))
        println("Node Links")
        println(node)
        println(dict_links_into_node)
        println(links_into_node)
        link_in[node] = links_into_node
        println(link_in[node])
        # Links emanating from a node
        dict_links_from_node = filter(kv -> (kv[2].from == node), network.links)
        links_from_node = sort(unique(keys(dict_links_from_node)))
        println(dict_links_from_node)
        println(links_from_node)
        link_out[node] = links_from_node
        println(link_out[node])
        # Links conected to a node (going to a node or emanating from a node)
        dict_links_tofrom_node = filter(kv -> ((kv[2].to == node) || (kv[2].from == node)), network.links)
        links_tofrom_node = sort(unique(keys(dict_links_tofrom_node)))
        link_tot[node] = links_tofrom_node
        println(links_tofrom_node)
        println(link_tot[node])
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

    # Maximum number of links connected to any node
    link_conn_max_length = maximum(unique(map(length,collect(values(link_tot)))))
    println("link_conn_max_length")
    println(link_conn_max_length)

    # Create a dictionary that maps each node to corresponding zones and the fraction of the node in each zone
    for node in node_set
        dict_node_zones = filter(kv ->(kv[2].location == node), network.zones)
        node_zones = [(x.zone, x.fraction) for x in values(dict_node_zones)]
        println("Node Zones")
        println(node)
        println(node_zones)
        dict_zones[node] = node_zones
        println(dict_zones[node])
    end
    # Create year, period and (year, period) sets for nodes with demands and existings
    demands_year_set = Vector{Year}()
    demands_period_set = Vector{PeriodName}()
    demands_year_period_set = Vector{Tuple{Year,PeriodName}}()
    existings_year_set = Vector{Year}()
    existings_period_set = Vector{PeriodName}()
    existings_year_period_set = Vector{Tuple{Year,PeriodName}}()

    for key in keys(demands)
	push!(demands_year_set, key.year)
	push!(demands_period_set, key.period)
        push!(demands_year_period_set, (key.year, key.period))
    end

    println("Existings years and periods")
    for value in values(existings)
        println(value.year, value.period)
    end

    for value in values(existings)
	for year in demands_year_set
	    if (year <= value.year)
               push!(existings_year_period_set, (value.year, value.period))
	    end
	end
    end


    demands_year_set = sort(unique(demands_year_set))
    demands_period_set = sort(unique(demands_period_set))
    demands_year_period_set = sort(unique(demands_year_period_set))
    existings_year_period_set = sort(unique(existings_year_period_set))

    println(" Demands Year Set")
    println(demands_year_set)
    println("Demands Period Set")
    println(demands_period_set)
    println("Demands Year Period Set")
    println(demands_year_period_set)
    println("Existings Year Period Set")
    println(existings_year_period_set)

    year_set = demands_year_set
    period_set = demands_period_set
    year_period_set = sort(unique(vcat(demands_year_period_set, existings_year_period_set)))

    # Create (year, period) set for nodes with storage
    storage_year_period_set = year_period_set               ####!!! Deepcopy?

    # Highest year
    max_year = maximum(map(x->x[1], year_period_set))
    # Final year = Highest year + 1
    fin_year = max_year + 1

    # Add (year, period) elements for final year to set of year-period tuples for storage
    for period in period_set
        push!(storage_year_period_set, (fin_year, period))
    end

    println("Year Period Set")
    println(year_period_set)

    println("Storage Year Period Set")
    println(storage_year_period_set)

    # Define dictionaries to be later used for defining JuMP cobstraints
    demands_year_dict = Dict{}()
    demands_period_dict = Dict{}()
    demands_year_period_dict = Dict{}()
    first_year_period_dict = Dict{}()
    int_year_period_dict = Dict{}()
    last_year_period_dict = Dict{}()
    allbutlast_year_period_dict = Dict{}()
    existings_year_dict = Dict{}()
    existings_period_dict = Dict{}()
    existings_year_period_dict = Dict{}()

    for node in demands_node_set
	node_keys = filter(k ->(k.location == node), keys(demands))
	node_year_set = sort(unique(getfield.(node_keys, :year)))
	node_period_set = sort(unique(getfield.(node_keys, :period)))
        node_year_period_set = Vector{Tuple{Year,PeriodName}}()
	for year in node_year_set
	    for period in node_period_set
                push!(node_year_period_set, (year, period))
	    end
	end
	demands_year_dict[node] = node_year_set
	demands_period_dict[node] = node_period_set
	demands_year_period_dict[node] = node_year_period_set
        if (length(demands_year_period_dict[node]) > 0)
	   first_year_period_dict[node] = [node_year_period_set[1]]
	   last_year_period_dict[node] = [node_year_period_set[end]]
	   allbutlast_year_period_dict[node] = node_year_period_set[1:end-1]
	   int_year_period_dict[node] = node_year_period_set[2:end-1]
	else
	   first_year_period_dict[node] = []
	   last_year_period_dict[node] = []
	   allbutlast_year_period_dict[node] = []
	   int_year_period_dict[node] = []
	end
    end

    for node in existings_node_set
	select_existings = filter(kv -> (kv[2].location == node), existings)
	node_year_set = sort(unique(getfield.(values(select_existings), :year)))
	node_period_set = sort(unique(getfield.(values(select_existings), :period)))
        node_year_period_set = Vector{Tuple{Year,PeriodName}}()
	for year in node_year_set
	    for period in node_period_set
                push!(node_year_period_set, (year, period))
	    end
	end
	existings_year_dict[node] = node_year_set
	existings_period_dict[node] = node_period_set
	existings_year_period_dict[node] = node_year_period_set
    end

    println("Demands_year_period_Dict")
    println(demands_year_period_dict)

    println("Existings_year_period_Dict")
    println(existings_year_period_dict)

    println("Year_period_dict")
    year_period_dict = merge(demands_year_period_dict, existings_year_period_dict)
    println(year_period_dict)

    println("Year_dict")
    year_dict = merge(demands_year_dict, existings_year_dict)
    println(year_dict)

    println("Period_dict")
    period_dict = merge(demands_period_dict, existings_period_dict)
    println(period_dict)

    # Reset node set to only those nodes with demands and existings
    if (length(unique(vcat(demands_node_set, existings_node_set))) < length(node_set))
       node_set = unique(vcat(demands_node_set, existings_node_set))
    end

    println("Modified Node Set")
    println(node_set)

    demands_total = Dict{Tuple{NetworkID,Tuple{Year,PeriodName}},Float64}()
    capacity = Dict{Tuple{NetworkID,Tuple{Year,PeriodName}},Float64}()
    existings_yield = Dict{Tuple{NetworkID,Technology,Tuple{Year,PeriodName}},Float64}()
    existings_cost = Dict{Tuple{NetworkID,Technology,Tuple{Year,PeriodName}},Float64}()
    pathway_technology = Dict{Tuple{Pathway,Int8},Technology}()
    technology_storage = Dict{Technology, Union{Technology,Nothing}}()

    for (key,value) in demands
	      demands_total[(key.location, (key.year, key.period))] = (value.fuelConsumption + value.nonfuelConsumption) * periods[key.period]
    end

    for value in values(existings)
	capacity[(value.location, (value.year, value.period))] = value.capacity * periods[value.period]
	existings_yield[(value.location, value.technology, (value.year, value.period))] = value.yield
	existings_cost[(value.location, value.technology, (value.year, value.period))] = value.cost
    end

    println("Existings Cost Dictionary")
    for kv in existings_cost
        println(kv)
    end

    for (key,value) in processlibrary.pathways
	pathway_technology[(key.pathway, key.stage)] = value.technology
    end

    println("Pathway Technology Dictionary")
    for kv in pathway_technology
        println(kv)
    end
    for (key,value) in processlibrary.costs
	      technology_storage[(key.technology)] = value.storage
    end

    println("Technology Storage Dictionary")
    for kv in technology_storage
        println(kv)
    end

    println("Total Demands")
    for (key,value) in demands_total
        println(key, demands_total[key])
    end


    # Define JuMP variables and constraints
    # JuMP production variable
    @variable(model, production[production_set, node_set, year_period_set] >= 0.0)
    # JuMP existing variable
    @variable(model, existingsvar[existings_set, existings_node_set, existings_year_period_set] >= 0.0)
    # JuMP consumption variable
    @variable(model, consumption[demands_node_set, demands_year_period_set] >= 0.0)

    # COMMENT: Why not consumption[j,t] <= demands_total[(j,t)]
    # JuMP constraint that consumption equals total demands
    @constraint(model, consumption_constr[j=demands_node_set, t=demands_year_period_dict[j]], consumption[j,t] == demands_total[(j,t)])

    # Capacity is not directly from processlibrary
    #@constraint(model, production_constr[i=production_set, j=node_set, t=year_period_dict[j]], production[i,j,t]<= capacity[(j,t)]) # Initially not implemented
    # COMMENT: Is the following line correct, because capacity is defined only for existings
    # JuMP constraint providing upper bound to production
	  @constraint(model, production_constr[i=production_set, j=existings_node_set, t=existings_year_period_dict[j]], production[i,j,t]<= capacity[(j,t)]) # Initially not implemented
    # COMMENT: Is existingsvar[i,j,t] == capacity[(j,t)]
    # JuMP constraint providing upper bound to existing
    @constraint(model, existings_constr[i=existings_set, j=existings_node_set, t=existings_year_period_dict[j]], existingsvar[i,j,t]<= capacity[(j,t)])

    # JuMP aggregate production variable
    @variable(model, aggregate_production[node_set, year_period_set] >= 0.0)

    # JuMP constraint that aggregate production equals sum of production and existing
    @constraint(model, aggregate_production_constr[j=node_set, t=year_period_dict[j]], aggregate_production[j,t] == sum(production[i,j,t] for i in production_set)
		+ sum(existingsvar[i,j,t] for i in existings_set if ((j in existings_node_set) && (t in existings_year_period_dict[j]))))

    # Unique pathway set
    pathway_set = sort(unique(getfield.(keys(processlibrary.pathways), :pathway)))
    println("Pathway Set")
    println(pathway_set)

    pathway_stages = Dict{}()
    pathways_extended = Dict{}()
    pathways_unextended = Dict{}()
    pathways_int_unextended = Dict{}()
    pathways_unextended_bfr = Dict{}()
    pathways_unextended_aft = Dict{}()
    path_unext_stage_w_storage = Dict{}()
    path_unext_stage_wo_storage = Dict{}()

    # Apportioning stages for each pathway into unextended, extended and other more detailed categories
    for pathway in pathway_set
	      pathway_keys = filter(k ->(k.pathway == pathway), keys(processlibrary.pathways))
        stage_set = sort(unique(getfield.(pathway_keys, :stage)))
	      println("Stage Set")
	      println(stage_set)

        pathway_stages[pathway] = maximum(stage_set)

	      dict_pathway_extended = filter(kv ->(kv[1].pathway == pathway && kv[2].extended == true), processlibrary.pathways)
        dict_pathway_unextended = filter(kv ->(kv[1].pathway == pathway && kv[2].extended == false), processlibrary.pathways)

        # Extended pathway stages
        pathways_extended[pathway] = sort(unique(getfield.(keys(dict_pathway_extended), :stage)))
        # Unextended pathway stages
        pathways_unextended[pathway] = sort(unique(getfield.(keys(dict_pathway_unextended), :stage)))
        # Intermediate Unextended pathway stages, i.e. not the first and ultimate stages
        pathways_int_unextended[pathway] = filter(e -> !(e in [1,pathway_stages[pathway]]),pathways_unextended[pathway])
        # Intermediate unextended pathway stages that are before the first extended stage
        pathways_unextended_bfr[pathway] = filter(e -> (e < minimum(pathways_extended[pathway])),pathways_int_unextended[pathway])
        # Intermediate unextended pathway stages that are after the final extended stage
        pathways_unextended_aft[pathway] = filter(e -> (e > maximum(pathways_extended[pathway])),pathways_int_unextended[pathway])
        println("Pathway Extended and Non-extended Stages")
        println(pathway)
        println(pathways_extended[pathway])
        println(pathways_unextended[pathway])
        println(pathways_int_unextended[pathway])
        println(pathways_unextended_bfr[pathway])
        println(pathways_unextended_aft[pathway])
    end

    # Maximum number of unextended pathway stages over all pathways
    println("Maximum Size of pathways_unextended dict over all pathways")
    pathways_unextended_max_length = last(sort(unique(map(length, collect(values(pathways_unextended))))))
    println(pathways_unextended_max_length)

    # Maximum number of extended pathway stages over all pathways
    println("Maximum Size of pathways_extended dict over all pathways")
    pathways_extended_max_length = last(sort(unique(map(length, collect(values(pathways_extended))))))
    println(pathways_extended_max_length)

    for w=1:length(pathway_set)
      for s=1:length(pathways_unextended[pathway_set[w]])
	for n=1:length(demands_node_set)
	  for t=1:length(first_year_period_dict[demands_node_set[n]])
	      println("Next")
	      println(pathway_set[w])
	      println(pathways_unextended[pathway_set[w]][s])
	      println(demands_node_set[n])
	      println(first_year_period_dict[demands_node_set[n]])
              print(next(first_year_period_dict[demands_node_set[n]][t], period_dict[demands_node_set[n]]))
	  end
	end
      end
    end

    # JuMP variables for unextended pathway stages
    @variable(model, pathway_details[p=1:length(pathway_set), s=pathways_unextended[pathway_set[p]], n=node_set, t=year_period_set] >= 0.0)
    println("Pathway Details")
    println(pathway_details)
    # JuMP variables for extended pathway stages with spatial flow between nodes
    @variable(model, pathway_edge_details[p=1:length(pathway_set), s=pathways_extended[pathway_set[p]], n=node_set, l=link_tot[n], t=year_period_set] >= 0.0)
    println("Pathway Edge Details")
    println(pathway_edge_details)
    # JuMP variables for unextended pathway stages
    @variable(model, pathway_details_storage[p=1:length(pathway_set), s=pathways_unextended[pathway_set[p]], n=node_set, t=storage_year_period_set] >= 0.0)
    # JuMP constraint for initial pathway stage (always unextended ??)
    @constraint(model, constr_init_stage[n=node_set, t=year_period_dict[n]], aggregate_production[n,t] == sum(pathway_details[p,1,n,t] for p=1:length(pathway_set)))
    # JuMP constraint for final pathway stage (always unextended ??)
    @constraint(model, constr_fin_stage[n=demands_node_set, t=demands_year_period_dict[n]], consumption[n,t] == sum(pathway_details[p,pathway_stages[pathway_set[p]],n,t] for p=1:length(pathway_set)))
    # JuMP constraint for intermediate unextended pathway stage before first extended stage for first time point
    @constraint(model, constr_int_unext_bfr_first[p=1:length(pathway_set), s=pathways_unextended_bfr[pathway_set[p]], n=demands_node_set, t=first_year_period_dict[n]], pathway_details_storage[p,s-1,n,next(t, period_dict[n])] + pathway_details[p,s,n,t] == pathway_details[p,s-1,n,t])
    # JuMP constraint for intermediate unextended pathway stage before first extended stage for intermediate time points
    @constraint(model, constr_int_unext_bfr_int[p=1:length(pathway_set), s=pathways_unextended_bfr[pathway_set[p]], n=demands_node_set, t=int_year_period_dict[n]], pathway_details_storage[p,s-1,n,next(t, period_dict[n])] + pathway_details[p,s,n,t] == pathway_details[p,s-1,n,t] + pathway_details_storage[p,s-1,n,prev(t, period_dict[n])])
    # JuMP constraint for intermediate unextended pathway stage before first extended stage for final time point
    @constraint(model, constr_int_unext_bfr_last[p=1:length(pathway_set), s=pathways_unextended_bfr[pathway_set[p]], n=demands_node_set, t=last_year_period_dict[n]], pathway_details[p,s,n,t] == pathway_details[p,s-1,n,t] + pathway_details_storage[p,s-1,n,prev(t, period_dict[n])])
    # JuMP constraint for intermediate unextended pathway stage after final extended stage for first time point
    @constraint(model, constr_int_unext_aft_first[p=1:length(pathway_set), s=pathways_unextended_aft[pathway_set[p]], n=demands_node_set, t=first_year_period_dict[n]], pathway_details_storage[p,s+1,n,next(t, period_dict[n])] + pathway_details[p,s,n,t] == pathway_details[p,s+1,n,t])
    # JuMP constraint for intermediate unextended pathway stage after final extended stage for intermediate time points
    @constraint(model, constr_int_unext_aft_int[p=1:length(pathway_set), s=pathways_unextended_aft[pathway_set[p]], n=demands_node_set, t=int_year_period_dict[n]], pathway_details_storage[p,s+1,n,next(t, period_dict[n])] + pathway_details[p,s,n,t] == pathway_details[p,s+1,n,t] + pathway_details_storage[p,s+1,n,prev(t, period_dict[n])])
    # JuMP constraint for intermediate unextended pathway stage after final extended stage for last time point
    @constraint(model, constr_int_unext_aft_last[p=1:length(pathway_set), s=pathways_unextended_aft[pathway_set[p]], n=demands_node_set, t=last_year_period_dict[n]], pathway_details[p,s,n,t] == pathway_details[p,s+1,n,t] + pathway_details_storage[p,s+1,n,prev(t, period_dict[n])])
    # JuMP constraint for extended pathway stages (always intermediate ??)
    #@constraint(model, constr_int_ext_first[p=1:length(pathway_set), s=pathways_extended[pathway_set[p]], n=node_set, t=first_year_period_dict[n]], pathway_details[p,s-1,n,t]
    @constraint(model, constr_int_ext[p=1:length(pathway_set), s=pathways_extended[pathway_set[p]], n=node_set, t=year_period_dict[n]], pathway_details[p,s-1,n,t]
			                                                                                                           + sum(pathway_edge_details[p,s,n,l,t] for l in link_in[n])
			                                                                                                           - sum(pathway_edge_details[p,s,other_node[(n,l)],l,t] for l in link_in[n] if other_node[(n,l)] in node_set)
			                                                                                                           + pathway_details_storage[p,s-1,n,t]
			                                                                                                           == pathway_details[p,s+1,n,t]
                                                                                                                                   #+ sum(pathway_edge_details[p,s,n,l,t] for l in link_out[n])
                                                                                                                                   #- sum(pathway_edge_details[p,s,other_node[(n,l)],l,t] for l in link_out[n]))
		                                                                                                                   - sum(pathway_edge_details[p,s,n,l,t] for l in link_out[n])
		                                                                                                                   + sum(pathway_edge_details[p,s,other_node[(n,l)],l,t] for l in link_out[n] if other_node[(n,l)] in node_set)
		                                                                                                                   + pathway_details_storage[p,s-1,n,next(t, period_dict[n])])

    #@constraint(model, constr_int_ext_last[p=1:length(pathway_set), s=pathways_extended[pathway_set[p]], n=node_set, t=last_year_period_dict[n]], pathway_details[p,s-1,n,t]
    #                                                                                                                                         + sum(pathway_edge_details[p,s,n,l,t] for l in link_in[n])
    #                                                                                                                                         - sum(pathway_edge_details[p,s,other_node[(n,l)],l,t] for l in link_in[n] if other_node[(n,l)] in node_set)
    #                                                                                                                                         + pathway_details_storage[p,s-1,n,t]
    #                                                                                                                                         == pathway_details[p,s+1,n,t]
    #                                                                                                                                         #+ sum(pathway_edge_details[p,s,n,l,t] for l in link_out[n])
    #                                                                                                                                         #- sum(pathway_edge_details[p,s,other_node[(n,l)],l,t] for l in link_out[n]))
    #                                                                                                                                         - sum(pathway_edge_details[p,s,n,l,t] for l in link_out[n])
    #                                                                                                                                         + sum(pathway_edge_details[p,s,other_node[(n,l)],l,t] for l in link_out[n] if other_node[(n,l)] in node_set))
    #@constraint(model, pathway_details_constr[i=1:length(pathway_set), j=1:pathway_stages[pathway_set[i]], k=node_set, t=year_period_set], pathway_details[i,j,k,t] <= stage_capacity[i,j,t]) #Initially not implemented
    #@constraint(model, pathway_details_constr[i=1:length(pathway_set), j=1:pathway_stages[pathway_set[i]], k=node_set, t=year_period_set], pathway_details[i,j,k,t] <= 50.0) #Initially not implemented

    production_construction = Array{Construction}(undef, niter, length(production_set), length(node_set), length(year_set))
    existings_existing = Array{Existing}(undef, niter, length(existings_set), length(existings_node_set), length(existings_year_period_set))
    pathway_details_construction = Array{Construction}(undef, niter, length(pathway_set), pathways_unextended_max_length, length(node_set), length(year_set))
    pathway_storage_construction = Array{Construction}(undef, niter, length(pathway_set), pathways_unextended_max_length, length(node_set), length(year_set))
    pathway_storage_dict = Array{Dict{}}(undef, niter, length(pathway_set), pathways_unextended_max_length, length(node_set), length(year_set))
    pathway_unext_input_dict = Array{Dict{}}(undef, niter, length(pathway_set), pathways_unextended_max_length, length(node_set), length(year_set))
    pathway_unext_output_dict = Array{Dict{}}(undef, niter, length(pathway_set), pathways_unextended_max_length, length(node_set), length(year_set))
    price_per_kg_of_hydrogen = Array{Float64}(undef, niter, length(pathway_set), pathways_unextended_max_length, length(node_set), length(year_set))
    pathway_edge_details_construction = Array{Construction}(undef, niter, length(pathway_set), pathways_extended_max_length, length(node_set), link_conn_max_length, length(year_set))
    pathway_edge_storage_construction = Array{Construction}(undef, niter, length(pathway_set), pathways_extended_max_length, length(node_set), link_conn_max_length, length(year_set))
    pathway_ext_input_dict = Array{Dict{}}(undef, niter, length(pathway_set), pathways_extended_max_length, length(node_set), link_conn_max_length, length(year_set))
    pathway_ext_output_dict = Array{Dict{}}(undef, niter, length(pathway_set), pathways_extended_max_length, length(node_set), link_conn_max_length, length(year_set))

    # Function to extract subset of data frame by variable, node, year and period
    function subdf0(DF,variable, node, year, period)
      sub0 = (flowDF.Variable .== variable) .&  (flowDF.Node .== node) .& (flowDF.Year .== year) .& (flowDF.Period .== period)
      return(sub0)
    end

    # Function to extract subset of data frame by variable, technology, node, and year
    function subdf1(DF,variable, technology, node, year)
      sub1 = (flowDF.Variable .== variable) .&  (flowDF.Technology .== technology) .& (flowDF.Node .== node) .& (flowDF.Year .== year)
      return(sub1)
    end

    # Function to extract subset of data frame by pathway and stage
    function subdf2(DF,pathway,stage)
      sub2 = (flowDF.Pathway .== pathway) .&  (flowDF.Stage .== stage)
      return(sub2)
    end

    # Function to extract subset of data frame by period and link
    function subdf3(DF,period,link)
      sub3 = (flowDF.Period .== period) .&  (flowDF.Link .== link)
      return(sub3)
    end

    # Initiate and populate data frame for flow variables before first iteration
    flowDF = DataFrame(Variable=String[], Technology=Technology[], Pathway=Pathway[] , Stage=String[] , Node=NetworkID[], Link=NetworkID[], Year=Year[], Period=PeriodName[], Flow=Float64[], Cost=Float64[])

    for p=1:length(production_set)
      for n=1:length(node_set)
        for t=1:length(year_dict[node_set[n]])
          for r=1:length(period_dict[node_set[n]])
            push!(flowDF, ("Production", production_set[p], "", "", node_set[n], "", year_dict[node_set[n]][t], period_dict[node_set[n]][r], 0.0, 0.0))
          end
        end
      end
    end

    for n in demands_node_set
      for t in demands_year_period_dict[n]
        push!(flowDF, ("Consumption", "", "", "", n, "", t[1], t[2], 0.0, 0.0))
      end
    end

    for p=1:length(existings_set)
      for n=1:length(existings_node_set)
        for t=1:length(existings_year_period_dict[existings_node_set[n]])
            push!(flowDF, ("Existing", existings_set[p], "", "", existings_node_set[n], "", existings_year_period_dict[existings_node_set[n]][t][1], existings_year_period_dict[existings_node_set[n]][t][2], capacity[existings_node_set[n], existings_year_period_dict[existings_node_set[n]][t]], existings_cost[existings_node_set[n], existings_set[p], existings_year_period_dict[existings_node_set[n]][t]] * capacity[existings_node_set[n], existings_year_period_dict[existings_node_set[n]][t]]))
        end
      end
    end

    for w=1:length(pathway_set)
      for s=1:length(pathways_unextended[pathway_set[w]])
        for n=1:length(node_set)
          for t=1:length(year_dict[node_set[n]])
            for r=1:length(period_dict[node_set[n]])
              if (node_set[n] in demands_node_set)
                 push!(flowDF, ("Pathway", pathway_technology[(pathway_set[w],pathways_unextended[pathway_set[w]][s])], pathway_set[w], string(pathways_unextended[pathway_set[w]][s]), node_set[n], "", year_dict[node_set[n]][t], period_dict[node_set[n]][r], demands_total[(node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r]))], 0.0))
              else
                 push!(flowDF, ("Pathway", pathway_technology[(pathway_set[w],pathways_unextended[pathway_set[w]][s])], pathway_set[w], string(pathways_unextended[pathway_set[w]][s]), node_set[n], "", year_dict[node_set[n]][t], period_dict[node_set[n]][r], 0.0, 0.0))
              end
            end
          end
        end
      end
    end

    for w=1:length(pathway_set)
      for s=1:length(pathways_extended[pathway_set[w]])
        for n=1:length(node_set)
          for l=1:length(link_tot[node_set[n]])
            for t=1:length(year_dict[node_set[n]])
              for r=1:length(period_dict[node_set[n]])
                if (node_set[n] in demands_node_set)
                   push!(flowDF, ("Pathway Edge", pathway_technology[(pathway_set[w],pathways_extended[pathway_set[w]][s])], pathway_set[w], string(pathways_extended[pathway_set[w]][s]), node_set[n], link_tot[node_set[n]][l], year_dict[node_set[n]][t], period_dict[node_set[n]][r], demands_total[(node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r]))], 0.0))
                else
                   push!(flowDF, ("Pathway Edge", pathway_technology[(pathway_set[w],pathways_extended[pathway_set[w]][s])], pathway_set[w], string(pathways_extended[pathway_set[w]][s]), node_set[n], link_tot[node_set[n]][l], year_dict[node_set[n]][t], period_dict[node_set[n]][r], 0.0, 0.0))
                end
              end
            end
          end
        end
      end
    end


    for w=1:length(pathway_set)
        path_unext_stage_wo_storage[pathway_set[w]] = []
        path_unext_stage_w_storage[pathway_set[w]] = []
        for s=1:length(pathways_unextended[pathway_set[w]]), n=1:length(node_set), t=1:length(year_dict[node_set[n]])
	    sub1 = subdf1(flowDF,"Pathway",pathway_technology[(pathway_set[w],pathways_unextended[pathway_set[w]][s])],node_set[n],year_dict[node_set[n]][t])
	    sub2 = subdf2(flowDF,pathway_set[w],string(pathways_unextended[pathway_set[w]][s]))
            capacity = maximum(flowDF[sub1 .& sub2, :Flow])
            distance = 0.0
            infrastructureid = string("I",'_',w,'_',s,'_',n,'_',t)
            dict_path_spec = filter(kv ->(kv[1].pathway == pathway_set[w] && kv[1].stage == pathways_unextended[pathway_set[w]][s]), processlibrary.pathways)
            technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
            storage = storagelookup(processlibrary, technology, year_dict[node_set[n]][t], capacity, distance, infrastructureid, node_set[n])

            if ((storage == "No Storage") || (storage == nothing))
               if !(pathways_unextended[pathway_set[w]][s] in path_unext_stage_wo_storage[pathway_set[w]])
                  push!(path_unext_stage_wo_storage[pathway_set[w]],pathways_unextended[pathway_set[w]][s])
               end
            else
               if !(pathways_unextended[pathway_set[w]][s] in path_unext_stage_w_storage[pathway_set[w]])
                  push!(path_unext_stage_w_storage[pathway_set[w]],pathways_unextended[pathway_set[w]][s])
               end
            end
        end
    end

    for w=1:length(pathway_set)
        println("Pathway Unextended Stages with Storage")
        println((pathway_set[w], path_unext_stage_w_storage[pathway_set[w]]))
        println("Pathway Unextended Stages without Storage")
        println((pathway_set[w], path_unext_stage_wo_storage[pathway_set[w]]))
    end

    pathway_set_w_node_storage = []
    for w=1:length(pathway_set)
        if (path_unext_stage_w_storage[pathway_set[w]] != [])
           push!(pathway_set_w_node_storage, pathway_set[w])
        end
    end

    println("Pathway Set with Node Storage")
    println(pathway_set_w_node_storage)

    w_map=Dict{}()
    s_map=Dict{}()

    for w=1:length(pathway_set_w_node_storage)
        w_map[w] = findfirst(isequal(pathway_set_w_node_storage[w]), pathway_set)
        for s=1:length(path_unext_stage_w_storage[pathway_set[w_map[w]]])
            s_map[(w_map[w],s)] = findfirst(isequal(path_unext_stage_w_storage[pathway_set[w_map[w]]][s]), pathways_unextended[pathway_set[w_map[w]]])
        end
    end


    for w=1:length(pathway_set_w_node_storage)
      for s=1:length(path_unext_stage_w_storage[pathway_set[w_map[w]]])
        for n=1:length(node_set)
          for t=1:length(year_dict[node_set[n]])
            for r=1:length(period_dict[node_set[n]])
                push!(flowDF, ("Storage", technology_storage[pathway_technology[(pathway_set[w_map[w]],path_unext_stage_w_storage[pathway_set[w_map[w]]][s])]], pathway_set[w_map[w]], string(path_unext_stage_w_storage[pathway_set[w_map[w]]][s]), node_set[n], "", year_dict[node_set[n]][t], period_dict[node_set[n]][r], 0.0, 0.0))
            end
          end
        end
      end
    end

    # JuMP constraint in which value of storage variable set to 0 for unextended pathway stages without storage
    @constraint(model, constr_unext_wo_storage[w=1:length(pathway_set), s=path_unext_stage_wo_storage[pathway_set[w]], n=node_set, t=year_period_dict[n]], pathway_details_storage[w,s,n,t] == 0.0)

    # JuMP constraint in which value of storage variable set to 0 for very first time point
    @constraint(model, constr_storage_first[w=1:length(pathway_set), s=path_unext_stage_w_storage[pathway_set[w]], n=node_set, t=first_year_period_dict[n]], pathway_details_storage[w,s,n,t] == 0.0)

    out_file = string("flowDF",0,".csv")
    CSV.write(out_file, flowDF)

    for i=1:niter
      for p=1:length(production_set), n=1:length(node_set), t=1:length(year_dict[node_set[n]])
          sub1 = subdf1(flowDF,"Production", production_set[p], node_set[n], year_dict[node_set[n]][t])
          capacity = maximum(flowDF[sub1, :Flow])
  	  distance = 0.0
          infrastructureid = string("I",'_',p,'_',n,'_',t)
          println("Iteration")
          println(i)
          production_construction[i,p,n,t] = processcosts(processlibrary, production_set[p], year_dict[node_set[n]][t], capacity, distance, infrastructureid, node_set[n])
      end

      for p=1:length(existings_set), n=1:length(existings_node_set), t=1:length(existings_year_period_dict[existings_node_set[n]])
          sub1 = subdf1(flowDF,"Existing",existings_set[p],existings_node_set[n],existings_year_period_dict[existings_node_set[n]][t][1])
          capacity = flowDF[sub1 .& (flowDF.Period .== existings_year_period_dict[existings_node_set[n]][t][2]), :Flow][1]
          cost = existings_cost[existings_node_set[n], existings_set[p], existings_year_period_dict[existings_node_set[n]][t]]
          println("Iteration")
          println(i)
          println("Flow")
          println(capacity)
          println("Cost")
  	      distance = 0.0
          println(cost)
          infrastructureid = string("I",'_',p,'_',n,'_',t)
          location = existings_node_set[n]
          technology = existings_set[p]
          year = existings_year_period_dict[existings_node_set[n]][t][1]
          period = existings_year_period_dict[existings_node_set[n]][t][2]
          yield = existings_yield[(location, technology, (year, period))]
          existings_existing[i,p,n,t] = Existing(location, technology, year, period, capacity, yield, cost)

      end

      for w=1:length(pathway_set), s=1:length(pathways_unextended[pathway_set[w]]), n=1:length(node_set), t=1:length(year_dict[node_set[n]])
          sub1 = subdf1(flowDF,"Pathway",pathway_technology[(pathway_set[w],pathways_unextended[pathway_set[w]][s])],node_set[n],year_dict[node_set[n]][t])
          sub2 = subdf2(flowDF,pathway_set[w],string(pathways_unextended[pathway_set[w]][s]))
          capacity = maximum(flowDF[sub1 .& sub2, :Flow])/maximum(sort(unique(values(periods))))
          println("Capacity")
          println(capacity)
  	      distance = 0.0
          infrastructureid = string("I",'_',w,'_',s,'_',n,'_',t)
          dict_path_spec = filter(kv ->(kv[1].pathway == pathway_set[w] && kv[1].stage == pathways_unextended[pathway_set[w]][s]), processlibrary.pathways)
          technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
          println("Technology before processcosts")
          println(technology)
          pathway_details_construction[i,w,s,n,t] = processcosts(processlibrary, technology, year_dict[node_set[n]][t], capacity, distance, infrastructureid, node_set[n])
          println("pathway_details_construction")
          println((pathway_set[w],pathways_unextended[pathway_set[w]][s],node_set[n],year_dict[node_set[n]][t], capacity, technology))
          println((i,w,s,n,t,pathway_details_construction[i,w,s,n,t]))
          #println(processinputs(processlibrary, technology, year_dict[node_set[n]][t], capacity, distance))
          #pathway_unext_input_dict[i,w,s,n,t] = processinputs(processlibrary, technology, year_dict[node_set[n]][t], capacity, distance)
          #pathway_unext_output_dict[i,w,s,n,t] = processoutputs(processlibrary, technology, year_dict[node_set[n]][t], capacity, distance)
          #println("Unextended Pathway Input Dictionary")
          #println((i,w,s,n,t,pathway_unext_input_dict[i,w,s,n,t]))
          #println("Unextended Pathway Output Dictionary")
          #println((i,w,s,n,t,pathway_unext_output_dict[i,w,s,n,t]))
          storage = storagelookup(processlibrary, technology, year_dict[node_set[n]][t], capacity, distance, infrastructureid, node_set[n])
          println("Storage")
          println(storage)
          if ((storage != "No Storage") && (storage != nothing))
             pathway_storage_construction[i,w,s,n,t] = processcosts(processlibrary, storage, year_dict[node_set[n]][t], capacity, distance, infrastructureid, node_set[n])
             pathway_storage_dict[i,w,s,n,t] = processinputs(processlibrary, storage, year_dict[node_set[n]][t], capacity, distance)
             println("Process Inputs for Storage")
             println(processinputs(processlibrary, storage, year_dict[node_set[n]][t], capacity, distance))
             materials = keys(pathway_storage_dict[i,w,s,n,t])
             println("Materials for Storage")
             println(materials)
             price_per_kg_of_hydrogen[i,w,s,n,t] = 0.0
             for material in materials
               # For this particular material, look up consumption/(kg of H2)
               quantity = pathway_storage_dict[i,w,s,n,t][material] * capacity
               price_material = 0.0
               println("Material")
               println(material)
               println("Material consumed per kg of H2")
               println(pathway_storage_dict[i,w,s,n,t][material])
               println("Total Material consumed")
               println(quantity)
               for k in 1:length(dict_zones[node_set[n]])
                   zone = dict_zones[node_set[n]][k][1]
                   fraction = dict_zones[node_set[n]][k][2]
                   println("Pricelookup Results")
                   println(zone)
                   println(fraction)
                   println(pricelookup(prices,usage,material,zone,year_dict[node_set[n]][t],quantity,true,true))
                   price_material = price_material + fraction * pricelookup(prices,usage,material,zone,year_dict[node_set[n]][t],quantity,true,true)
               end
               println("Material Price in \$/unit")
               println(price_material,",",quantity,",",capacity)
               price_per_kg_of_hydrogen[i,w,s,n,t] = price_per_kg_of_hydrogen[i,w,s,n,t]  + price_material * quantity/capacity
             end
             println("price per kg of hydrogen")
             println(price_per_kg_of_hydrogen[i,w,s,n,t])
             println("pathway_storage_construction")
             println((i,w,s,n,t,pathway_storage_construction[i,w,s,n,t]))
          end
      end

      for w=1:length(pathway_set), s=1:length(pathways_extended[pathway_set[w]]), n=1:length(node_set), l=1:length(link_tot[node_set[n]]), t=1:length(year_dict[node_set[n]])
          sub1 = subdf1(flowDF, "Pathway Edge", pathway_technology[(pathway_set[w],pathways_extended[pathway_set[w]][s])], node_set[n], year_dict[node_set[n]][t])
          sub2 = subdf2(flowDF, pathway_set[w], string(pathways_extended[pathway_set[w]][s]))
          sub3 = (flowDF.Link .== link_tot[node_set[n]][l])
          capacity = maximum(flowDF[sub1 .& sub2 .& sub3, :Flow])
  	      distance = link_length[link_tot[node_set[n]][l]]
          link_key = link_tot[node_set[n]][l]
          println("Iteration")
          println(i)
          println("Node")
          println(node_set[n])
          println("To/From")
          println(network.links[link_key].from)
          println(network.links[link_key].to)
          println("Distance")
          println(distance)
          println("Capacity")
          println(capacity)
          println("Pathway Technology")
          println(pathway_technology[(pathway_set[w],pathways_extended[pathway_set[w]][s])])

          infrastructureid = string("I",'_',w,'_',n,'_',l,'_',t)
          dict_path_spec = filter(kv ->(kv[1].pathway == pathway_set[w] && kv[1].stage == pathways_extended[pathway_set[w]][s]), processlibrary.pathways)
          technology = sort(unique(getfield.(values(dict_path_spec), :technology)))[1]
          println("Technology")
          println(technology)
          pathway_edge_details_construction[i,w,s,n,l,t] = processcosts(processlibrary, technology, year_dict[node_set[n]][t], capacity, distance, infrastructureid, node_set[n])
          println("pathway_edge_details_construction")
          println(pathway_edge_details_construction[i,w,s,n,l,t])
      end


      println("List of pathway storage constructions")
      for w=1:length(pathway_set_w_node_storage), s=1:length(path_unext_stage_w_storage[pathway_set[w_map[w]]]), n=1:length(node_set), t=1:length(year_dict[node_set[n]]), r=1:length(period_dict[node_set[n]])
          println((w,s,n,t,r, pathway_storage_construction[i,w_map[w],s_map[(w_map[w],s)],n,t], pathway_details_storage[w_map[w], path_unext_stage_w_storage[pathway_set[w_map[w]]][s], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])]), price_per_kg_of_hydrogen[i,w_map[w],s_map[(w_map[w],s)],n,t])
      end

      # Define objective function and solve optimization model

      println("Solve model")
      o=@objective(model,Min,sum(production_construction[i,p,n,t].variablecost * production[production_set[p], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])]
                             for p=1:length(production_set), n=1:length(node_set), t=1:length(year_dict[node_set[n]]), r=1:length(period_dict[node_set[n]])) +
                             sum(existings_existing[i,p,n,t].cost * existingsvar[existings_set[p], existings_node_set[n], existings_year_period_dict[existings_node_set[n]][t]]
                             for p=1:length(existings_set), n=1:length(existings_node_set), t=1:length(existings_year_period_dict[existings_node_set[n]])) +
                             sum(pathway_details_construction[i,w,s,n,t].variablecost * pathway_details[w, pathways_unextended[pathway_set[w]][s], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])]
                             for w=1:length(pathway_set), s=1:length(pathways_unextended[pathway_set[w]]), n=1:length(node_set), t=1:length(year_dict[node_set[n]]), r=1:length(period_dict[node_set[n]])) +
			                       sum((pathway_storage_construction[i,w_map[w],s_map[(w_map[w],s)],n,t].variablecost + price_per_kg_of_hydrogen[i,w_map[w],s_map[(w_map[w],s)],n,t]) * pathway_details_storage[w_map[w], path_unext_stage_w_storage[pathway_set[w_map[w]]][s], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])]
                             for w=1:length(pathway_set_w_node_storage), s=1:length(path_unext_stage_w_storage[pathway_set[w]]), n=1:length(node_set), t=1:length(year_dict[node_set[n]]), r=1:length(period_dict[node_set[n]])) +
                             sum(pathway_edge_details_construction[i,w,s,n,l,t].variablecost * pathway_edge_details[w, pathways_extended[pathway_set[w]][s], node_set[n], link_tot[node_set[n]][l], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])]
                             for w=1:length(pathway_set), s=1:length(pathways_extended[pathway_set[w]]), n=1:length(node_set), l=1:length(link_tot[node_set[n]]), t=1:length(year_dict[node_set[n]]), r=1:length(period_dict[node_set[n]]))
      )
      println("Objective Function Start")
      println(o)
      println("Objective Function End")
      optimize!(model)
      println(termination_status(model))

      # Update data frame for each optimization iteration

      for p=1:length(production_set)
        for n=1:length(node_set)
          for t=1:length(year_dict[node_set[n]])
            for r=1:length(period_dict[node_set[n]])
              println("Iteration ",i," Optimized Production")
              sub1 = subdf1(flowDF, "Production",production_set[p],node_set[n],year_dict[node_set[n]][t])
              sub2 = (flowDF.Period .== period_dict[node_set[n]][r])
              flowDF[sub1 .& sub2, :Flow] .= value(production[production_set[p], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])])
              flowDF[sub1 .& sub2, :Cost] .= production_construction[i,p,n,t].variablecost * flowDF[sub1 .& sub2, :Flow]
            end
          end
        end
      end

      println(flowDF)
      for n in demands_node_set
        for t in demands_year_period_dict[n]
          sub0 = subdf0(flowDF,"Consumption",n,t[1],t[2])
          flowDF[sub0, :Flow] .= value(consumption[n,t])
        end
      end


      for p=1:length(existings_set)
        for n=1:length(existings_node_set)
          for t=1:length(existings_year_period_dict[existings_node_set[n]])
                println("Iteration ",i)
                println("Optimized Existing")
                sub1 = subdf1(flowDF,"Existing",existings_set[p],existings_node_set[n],existings_year_period_dict[existings_node_set[n]][t][1])
                sub2 = (flowDF.Period .== existings_year_period_dict[existings_node_set[n]][t][2])
                flowDF[sub1 .& sub2, :Flow] .= value(existingsvar[existings_set[p], existings_node_set[n], existings_year_period_dict[existings_node_set[n]][t]])
	        flowDF[sub1 .& sub2, :Cost] .= existings_existing[i,p,n,t].cost * flowDF[sub1 .& sub2, :Flow]
          end
        end
      end

      for w=1:length(pathway_set)
        for s=1:length(pathways_unextended[pathway_set[w]])
          for n=1:length(node_set)
            for t=1:length(year_dict[node_set[n]])
              for r=1:length(period_dict[node_set[n]])
                sub1 = subdf1(flowDF,"Pathway",pathway_technology[(pathway_set[w],pathways_unextended[pathway_set[w]][s])],node_set[n],year_dict[node_set[n]][t])
                sub2 = subdf2(flowDF,pathway_set[w],string(pathways_unextended[pathway_set[w]][s]))
                sub3 = (flowDF.Period .== period_dict[node_set[n]][r])
                flowDF[sub1 .& sub2 .& sub3, :Flow] .= value(pathway_details[w, pathways_unextended[pathway_set[w]][s], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])])
	        flowDF[sub1 .& sub2 .& sub3, :Cost] .= pathway_details_construction[i,w,s,n,t].variablecost * flowDF[sub1 .& sub2 .& sub3, :Flow]
              end
            end
          end
        end
      end


      for w=1:length(pathway_set)
        for s=1:length(pathways_extended[pathway_set[w]])
          for n=1:length(node_set)
            for l=1:length(link_tot[node_set[n]])
              for t=1:length(year_dict[node_set[n]])
                for r=1:length(period_dict[node_set[n]])
                  sub1 = subdf1(flowDF,"Pathway Edge",pathway_technology[(pathway_set[w],pathways_extended[pathway_set[w]][s])],node_set[n],year_dict[node_set[n]][t])
                  sub2 = subdf2(flowDF,pathway_set[w],string(pathways_extended[pathway_set[w]][s]))
                  sub3 = subdf3(flowDF,period_dict[node_set[n]][r],link_tot[node_set[n]][l])
                  flowDF[sub1 .& sub2 .& sub3, :Flow] .= value(pathway_edge_details[w, pathways_extended[pathway_set[w]][s], node_set[n], link_tot[node_set[n]][l], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])])
                  flowDF[sub1 .& sub2 .& sub3, :Cost] .= pathway_edge_details_construction[i,w,s,n,l,t].variablecost * flowDF[sub1 .& sub2 .& sub3, :Flow]
                end
              end
            end
          end
        end
      end

      for w=1:length(pathway_set_w_node_storage)
        for s=1:length(path_unext_stage_w_storage[pathway_set[w_map[w]]])
          for n=1:length(node_set)
            for t=1:length(year_dict[node_set[n]])
              for r=1:length(period_dict[node_set[n]])
                println("Storage Flow")
                println(value(pathway_details_storage[w_map[w], path_unext_stage_w_storage[pathway_set[w_map[w]]][s], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])]))
                sub1 = subdf1(flowDF,"Storage",pathway_storage_construction[i,w_map[w],s_map[(w_map[w],s)],n,t].technology,node_set[n],year_dict[node_set[n]][t])
                sub2 = subdf2(flowDF,pathway_set[w_map[w]],string(path_unext_stage_w_storage[pathway_set[w_map[w]]][s]))
                sub3 = (flowDF.Period .== period_dict[node_set[n]][r])
                flowDF[sub1 .& sub2 .& sub3, :Flow] .= abs(value(pathway_details_storage[w_map[w], path_unext_stage_w_storage[pathway_set[w_map[w]]][s], node_set[n], (year_dict[node_set[n]][t], period_dict[node_set[n]][r])]) - value(pathway_details_storage[w_map[w], path_unext_stage_w_storage[pathway_set[w_map[w]]][s], node_set[n], next((year_dict[node_set[n]][t], period_dict[node_set[n]][r]), period_dict[node_set[n]])]))
                flowDF[sub1 .& sub2 .& sub3, :Cost] .= (pathway_storage_construction[i,w_map[w],s_map[(w_map[w],s)],n,t].variablecost + price_per_kg_of_hydrogen[i,w_map[w],s_map[(w_map[w],s)],n,t]) * flowDF[sub1 .& sub2 .& sub3, :Flow]
              end
            end
          end
        end
      end

      print(flowDF)

      # Output data frame for each iteration

      out_file = string("flowDF",i,".csv")
      CSV.write(out_file, flowDF)
    end

    return(model)
end



end
