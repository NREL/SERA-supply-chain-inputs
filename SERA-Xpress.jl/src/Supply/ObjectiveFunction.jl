function create_objective_function(annualization,
                                   i,
                                   model,
                                   new_production_cap,pathway_unextended_cap, pathway_storage_cap, pathway_extended_cap,
                                   new_production, existing_production, pathway_stage_flow, pathway_stage_stor_in, pathway_stage_stor_out, link_flow,
                                   node_set, node_existings, production_nodes, link_tot, node_production_set, pathway_set,
                                   pathways_unextended, pathway_storage_stages,  pathways_extended,
                                   year_set, year_period_set, total_year_set, npv_array,
                                   yearly_cap_cost_production, yearly_cap_cost_unextended, yearly_cap_cost_storage, yearly_cap_cost_extended,
                                   yearly_fixed_cost_production, yearly_fixed_cost_unextended, yearly_fixed_cost_storage, yearly_fixed_cost_extended,
                                   production_data, production_price_per_kg_of_hydrogen,
                                   existings_data, existings_price_per_kg_of_hydrogen,
                                   pathway_unext_flow_data, pathway_unextended_price_per_kg_of_hydrogen,
                                   pathway_unext_storage_data, storage_price_per_kg_of_hydrogen,
                                   pathway_ext_link_data, pathway_extended_price_per_kg_of_hydrogen
                                   )

    if annualization
        o=@objective(model, Min,
        #Capital Costs
                                sum(npv_array[t] * (
                                    sum(yearly_cap_cost_production[n, x, y, t] for n in production_nodes, x in node_production_set[n]) +
                                    sum(yearly_cap_cost_unextended[n, p, s, y, t] for n in node_set, p in pathway_set, s in pathways_unextended[p]) +
                                    sum(yearly_cap_cost_storage[n, p, s, y, t] for n in node_set, p in pathway_set, s in pathway_storage_stages[p]) +
                                    sum(yearly_cap_cost_extended[n, l, p, s, y, t] for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p])
                                )
                                    for y in year_set, t in year_set)+

        #Fixed Costs
                                    sum(npv_array[t] * (
                                    sum(yearly_fixed_cost_production[n, x, y, t] for n in production_nodes, x in node_production_set[n]) +
                                    sum(yearly_fixed_cost_unextended[n, p, s, y, t] for n in node_set, p in pathway_set, s in pathways_unextended[p]) +
                                    sum(yearly_fixed_cost_storage[n, p, s, y, t] for n in node_set, p in pathway_set, s in pathway_storage_stages[p]) +
                                    sum(yearly_fixed_cost_extended[n, l, p, s, y, t] for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p])
                                                    )
                                for y in year_set, t in year_set) +
        #Operation Costs
                                sum(npv_array[t[1]] * (
                                    sum((production_data[i, n, x, string(t[1])].variablecost + production_price_per_kg_of_hydrogen[n, x, t]) * new_production[n, x, y, t] for n in production_nodes, x in node_production_set[n], y in total_year_set) +
                                    sum((existings_data[n, x, t].cost + existings_price_per_kg_of_hydrogen[n, x, t]) * existing_production[n, x, t] for n in node_set, x in node_existings[n]) +
                                    sum((pathway_unext_flow_data[i, n, p, s, string(t[1])].variablecost + pathway_unextended_price_per_kg_of_hydrogen[n, p, s, t]) * pathway_stage_flow[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathways_unextended[p], y in total_year_set) +
                                    sum((pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[n, p, s, t]) * pathway_stage_stor_in[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in total_year_set) +
                                    sum((pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[n, p, s, t]) * pathway_stage_stor_out[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in total_year_set) +
                                    sum((pathway_ext_link_data[i, n, p, s, l, string(t[1])].variablecost + pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, t]) * link_flow[n, l, p, s, y, t]
                                        for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in total_year_set)
                                                        )
                                for t in year_period_set))
    else
        o=@objective(model, Min,
        #Capital Costs
                                sum(npv_array[t] * (
                                    sum(new_production_cap[n, x, t] * production_data[i, n, x, string(t)].capitalcost  for n in production_nodes, x in node_production_set[n]) +
                                    sum(pathway_unextended_cap[n, p, s, t] * pathway_unext_flow_data[i, n, p, s, string(t)].capitalcost for n in node_set, p in pathway_set, s in pathways_unextended[p]) +
                                    sum(pathway_storage_cap[n, p, s, t] * pathway_unext_storage_data[i, n, p, s, string(t)].capitalcost for n in node_set, p in pathway_set, s in pathway_storage_stages[p]) +
                                    sum(pathway_extended_cap[n, l, p, s, t] * pathway_ext_link_data[i, n, p, s, l, string(t)].capitalcost for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p])
                                )
                                    for t in year_set)+

        #Fixed Costs
                                    sum(npv_array[t] * (
                                    sum(yearly_fixed_cost_production[n, x, y, t] for n in production_nodes, x in node_production_set[n]) +
                                    sum(yearly_fixed_cost_unextended[n, p, s, y, t] for n in node_set, p in pathway_set, s in pathways_unextended[p]) +
                                    sum(yearly_fixed_cost_storage[n, p, s, y, t] for n in node_set, p in pathway_set, s in pathway_storage_stages[p]) +
                                    sum(yearly_fixed_cost_extended[n, l, p, s, y, t] for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p])
                                                    )
                                for y in year_set, t in year_set) +
        #Operation Costs
                                sum(npv_array[t[1]] * (
                                    sum((production_data[i, n, x, string(t[1])].variablecost + production_price_per_kg_of_hydrogen[n, x, t]) * new_production[n, x, y, t] for n in production_nodes, x in node_production_set[n], y in total_year_set) +
                                    sum((existings_data[n, x, t].cost + existings_price_per_kg_of_hydrogen[n, x, t]) * existing_production[n, x, t] for n in node_set, x in node_existings[n]) +
                                    sum((pathway_unext_flow_data[i, n, p, s, string(t[1])].variablecost + pathway_unextended_price_per_kg_of_hydrogen[n, p, s, t]) * pathway_stage_flow[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathways_unextended[p], y in total_year_set) +
                                    sum((pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[n, p, s, t]) * pathway_stage_stor_in[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in total_year_set) +
                                    sum((pathway_unext_storage_data[i, n, p, s, string(t[1])].variablecost + storage_price_per_kg_of_hydrogen[n, p, s, t]) * pathway_stage_stor_out[n, p, s, y, t]
                                        for n in node_set, p in pathway_set, s in pathway_storage_stages[p], y in total_year_set) +
                                    sum((pathway_ext_link_data[i, n, p, s, l, string(t[1])].variablecost + pathway_extended_price_per_kg_of_hydrogen[n, p, s, l, t]) * link_flow[n, l, p, s, y, t]
                                        for n in node_set, l in link_tot[n], p in pathway_set, s in pathways_extended[p], y in total_year_set)
                                                        )
                                for t in year_period_set))
    end
    return o
end
