"""
Look up and interpolate prices.

FIXME: document arguments and add example.
"""

module LookupFunctions

using CSV
using DataFrames
using DataFramesMeta
using SERA.Supply.Types
using SERA.Util

import SERA.Supply.Types: MaterialName, MaterialUsage, Prices, ZoneID, Year

export pricelookup
export pricefunction
export processcosts
export processinputs
export processoutputs
export storagelookup
export extrap_cost
export find_init_capacity_and_utilization
# Function to find two values from array nearest to given value
function findtwonearest(x, xlist)
    x1 = xlist[findmin(abs.(xlist .- x))[2]]
    deleteat!(xlist, findmin(abs.(xlist .- x))[2])
    if length(xlist) > 0
       x2 = xlist[findmin(abs.(xlist .- x))[2]]
    else
       x2 = x1
    end
    return x1,x2
end

function findnearestleq(x, xlist)
    @assert x >= minimum(xlist)
    x1 = xlist[findmin(abs.(xlist .- x))[2]]
    if x1 <= x
       return x1
    else
       ind = findfirst(isequal(x1), xlist)
       return xlist[ind-1]
    end
end

# Function for linear interpolation
function linearinterpolation(x, x1, x2, y1, y2)
    y = y1 + (x - x1)*(y2 - y1)/(x2 - x1)
    return y
end


# Function to Estimate price for quantity for given material, zone and year
function pqinterpolate(prices, material, zone, year, period, quantity, pq_interpolate=true)
   priceqty = Dict{Float64,Float64}()
   for key in keys(prices[PriceKey(material,zone,year,period)])
      priceqty[key] = prices[PriceKey(material,zone,year,period)][key].price
   end
   qtylist = collect(keys(priceqty))
   qty1,qty2 = findtwonearest(quantity, qtylist)
   price1 = priceqty[qty1]
   if ((abs(qty1 - qty2) < 1.0e-8) | (qty1 == qty2)) # check if qty1 and qty2 are equal
      pq_interpolate = false
   end
   # Within a given year, decide if to find price for quantity by interpolation or choose price corresponding to nearest quantity
   if pq_interpolate
      price2 = priceqty[qty2]
      price = linearinterpolation(quantity, qty1, qty2, price1, price2)
   else
      price = price1
   end
   return price
end


function pricelookup(prices :: Prices, usage :: MaterialUsage, material :: MaterialName, zone :: ZoneID, year :: Year, period:: String, quantity :: Float64, pq_interpolate=false, py_interpolate=true) :: Float64
    # First,  find how much material is already used in the zone and year.

    if PriceKey(material,zone,year,period) in keys(usage)
      materialusage = usage[PriceKey(material,zone,year,period)]

      # Second, the new quantity used will be that plus the amount being priced.
      quantity = quantity + materialusage
      # Third, find the two points on the price curve that bracket the year.

      # Create list of all years for which price-quantity supply curves are available for given material and zone
      selected_prices = filter(kv->kv[1].material == material && kv[1].zone == zone, prices)
      years = sort(unique(getfield.(keys(selected_prices), :year)))
      periods = sort(unique(getfield.(keys(selected_prices), :period)))

      # Check if price-quantity supply curve is available for given material in a given zone in a given year
      if year in years && period in periods
         price = pqinterpolate(prices, material, zone, year, period, quantity, pq_interpolate)
      else
         #=
         # In case of no price-quantity supply curve for a given year, find the prices for quantity for two nearest years

         year1,year2 = findtwonearest(year, years)
         priceyear1 = pqinterpolate(prices, material, zone, year1, period, quantity, pq_interpolate)
         priceyear2 = pqinterpolate(prices, material, zone, year2, period, quantity, pq_interpolate)
         # Decide if to find price by interpolation from prices for two nearest years or just choose price for nearest year
         if py_interpolate
            price = linearinterpolation(year, year1, year2, price1, price2)
         else
            priceyear = Dict()
            priceyear[year1] = price1
            priceyear[year2] = price2
            nearest_year = findtwonearest(year, collect(keys(priceyear)))[1]
            price = priceyear[nearest_year]
         end
         =#

         #In case of no price-quantity supply curve for a given year and period, give an error
         error("Price for material $(material) not for $(period), $(year)")
      end
      return price
   else
      return 0.0
   end
end

"""
Return a function for looking up prices.

FIXME: document arguments and add example.
"""
function pricefunction(prices :: Prices, usage :: MaterialUsage, material :: MaterialName, zone :: ZoneID, pq_interpolate=true, py_interpolate=true)
    function f(year :: Year, quantity :: Float64)
        pricelookup(prices, usage, material, zone, pq_interpolate, py_interpolate)
    end
    return
end

function find_init_capacity_and_utilization(library :: ProcessLibrary, technology :: Technology, year :: Year)

   selected_costs = filter(kv->kv[1].technology == technology, library.costs)

   years = sort(unique(getfield.(keys(selected_costs), :year)))

   year_less= findnearestleq(year, years)

   selected_year_costs = filter(kv->kv[1].year == year_less, selected_costs)

    num_entries = length(keys(selected_year_costs))

    if num_entries == 0
      capacity = 0.0
      utilization = 0.0
    elseif num_entries <= 2
      capacity = sort(unique(getfield.(keys(selected_year_costs), :nameplate)))[1]
      utilization = sort(unique(getfield.(keys(selected_year_costs), :dutycycle)))[1]
    else
      capacity = sort(unique(getfield.(keys(selected_year_costs), :nameplate)))[end]
      utilization = sort(unique(getfield.(keys(selected_year_costs), :dutycycle)))[1]
    end

    return capacity, utilization
end

function processcosts(library :: ProcessLibrary, technology :: Technology, year :: Year, capacity :: Float64, distance :: Float64, infrastructureid :: InfrastructureID, networkid :: NetworkID) :: Construction
    #=
    First find the nearest records in `library.costs`:
    The technology name must match exactly.
    =#
    selected_costs = filter(kv->kv[1].technology == technology, library.costs)

    # Use the greatest year in the process keys that is less than or equal to `year`.
    years = sort(unique(getfield.(keys(selected_costs), :year)))

    year_less= findnearestleq(year, years)

    selected_year_costs = filter(kv->kv[1].year == year_less, selected_costs)

    #=
    Compare `nameplate * dutycycle` to `capacity`:
    if `capacity` is smaller than the smallest `nameplate * dutycycle`, then just use that key.
    otherwise, find the key with the greatest `nameplate * dutycycle` that is less than or equal to `capacity`.
    Let `key` be the key and `value` be the value.
    if capacity less than the smallest `nameplate * dutycycle`
    then ratio = 1
    otherwise ratio = capacity / key.nameplate / key.dutycycle
    =#

    dict_effcap = Dict()
    for key in keys(selected_year_costs)
        effective_capacity = key.nameplate * key.dutycycle
        dict_effcap[effective_capacity] = key
    end
    list_effcap = sort(collect(keys(dict_effcap)))

    lowest_effcap = list_effcap[1]

    if (capacity < lowest_effcap)
       selected_key = dict_effcap[lowest_effcap]
       ratio = 1
    else
       effcap_less = findnearestleq(capacity, list_effcap)
       selected_key = dict_effcap[effcap_less]
       ratio = capacity / selected_key.nameplate / selected_key.dutycycle
    end
    selected_value = selected_year_costs[selected_key]

    #The return value should be Construction structure with the following:

    #result.infrastructureid = infrastructureid
    #result.networkid = networkid
    technology = selected_key.technology
    productive = selected_value.productive
    #result.year = year
    lifetime = selected_value.lifetime
    nameplate = selected_key.nameplate / selected_key.dutycycle
    dutycycle = selected_key.dutycycle
    #result.distance = distance
    cap_cost_polynomial = selected_value.capitalcost
    base_cap_cost = cap_cost_polynomial.constant +
                      (cap_cost_polynomial.linear * distance) +
                      (cap_cost_polynomial.quadratic * (distance ^ 2))

    capitalcost = base_cap_cost * ratio^selected_value.scaling / (nameplate * ratio)

    base_fixed_cost = base_cap_cost * selected_value.fixedcostfrac
    fixedcost = base_fixed_cost * ratio^selected_value.scaling / (nameplate * ratio)

    variablecost = selected_value.variablecost + distance * selected_value.variablecoststretch

    storage_polynomial = selected_value.storage.capacitypolynomial
    storagecapacity = storage_polynomial.constant +
                      (storage_polynomial.linear * distance) +
                      (storage_polynomial.quadratic * (distance ^ 2)) +
                      (storage_polynomial.cubic * (distance ^ 3))

    result = Construction(infrastructureid, networkid, technology, productive, year, lifetime, nameplate,
                          dutycycle, distance, capitalcost, fixedcost, variablecost, storagecapacity)

    return result

end


function processinputs(library :: ProcessLibrary, technology :: Technology, year :: Year, capacity :: Float64, distance :: Float64) :: Dict{MaterialName,Float64}
    #=
    The same lookup occurs as in `processcosts`, except dutycycle is ignored.
    The output is a dictionary for where the key is each material and the value is `consumptionrate + distance * consumptionratestretch`

    First find the nearest records in `library.inputs`:
    The technology name must match exactly.
    =#

    selected_inputs = filter(kv->kv[1].technology == technology, library.inputs)

    if isempty(selected_inputs)
      return Dict{MaterialName,Float64}()
    else
      # Use the greatest year in the process keys that is less than or equal to `year`.
      years = sort(unique(getfield.(keys(selected_inputs), :year)))
      year_less= findnearestleq(year, years)
      selected_year_inputs = filter(kv->kv[1].year == year_less, selected_inputs)
      dict_effcap = Dict()
      for key in keys(selected_year_inputs)
         effective_capacity = key.nameplate
         dict_effcap[key] = effective_capacity
      end
      list_effcap = sort(collect(values(dict_effcap)))
      lowest_effcap = list_effcap[1]
      if (capacity < lowest_effcap)
         selected_dict_effcap = filter(kv->kv[2] == lowest_effcap, dict_effcap)
      else
         effcap_less = findnearestleq(capacity, list_effcap)
         selected_dict_effcap = filter(kv->kv[2] == effcap_less, dict_effcap)
      end

      dict = Dict{MaterialName,Float64}()
      for key in keys(selected_dict_effcap)
         value = selected_year_inputs[key]
         dict[key.material] = value.consumptionrate + distance * value.consumptionratestretch
      end

      return dict
   end

end


function processoutputs(library :: ProcessLibrary, technology :: Technology, year :: Year, capacity :: Float64, distance :: Float64) :: Dict{MaterialName,Float64}
    #=
    The same lookup occurs as in `processcosts`, except dutycycle is ignored.
    the output is a dictionary for where the key is each material and the value is `productionrate + distance * productionratestretch`

    First find the nearest records in `library.outputs`:
    The technology name must match exactly.
    =#
    selected_outputs = filter(kv->kv[1].technology == technology, library.outputs)
    if isempty(selected_outputs)
      return Dict{MaterialName,Float64}()
    else
      # Use the greatest year in the process keys that is less than or equal to `year`.
      years = sort(unique(getfield.(keys(selected_outputs), :year)))
      year_less= findnearestleq(year, years)
      selected_year_outputs = filter(kv->kv[1].year == year_less, selected_outputs)

      dict_effcap = Dict()
      for key in keys(selected_year_outputs)
         effective_capacity = key.nameplate
         dict_effcap[key] = effective_capacity
      end
      #list_effcap = sort(collect(values(dict_effcap)))
      list_effcap = sort(unique(collect(values(dict_effcap))))
      lowest_effcap = list_effcap[1]
      if (capacity < lowest_effcap)
         selected_diff_effcap = filter(kv-> kv[2] == lowest_effcap, dict_effcap)
      else
         effcap_less = findnearestleq(capacity, list_effcap)
         selected_diff_effcap = filter(kv-> kv[2] == effcap_less, dict_effcap)
      end

      dict = Dict{MaterialName,Float64}()
      for key in keys(selected_diff_effcap)
         value = selected_year_outputs[key]
         dict[key.material] = value.productionrate + distance * value.productionratestretch
      end

      return dict
   end

end

function storagelookup(library :: ProcessLibrary, technology :: Technology, year :: Year, capacity :: Float64, distance :: Float64, infrastructureid :: InfrastructureID, networkid :: NetworkID) :: Union{Storage,Nothing}
    #=
    First find the nearest records in `library.costs`:
    The technology name must match exactly.
    =#
    #println(library.costs)
    selected_costs = filter(kv->kv[1].technology == technology, library.costs)

    # Use the greatest year in the process keys that is less than or equal to `year`.
    years = sort(unique(getfield.(keys(selected_costs), :year)))
    year_less= findnearestleq(year, years)
    selected_year_costs = filter(kv->kv[1].year == year_less, selected_costs)

    #=
    Compare `nameplate * dutycycle` to `capacity`:
    if `capacity` is smaller than the smallest `nameplate * dutycycle`, then just use that key.
    otherwise, find the key with the greatest `nameplate * dutycycle` that is less than or equal to `capacity`.
    Let `key` be the key and `value` be the value.
    if capacity less than the smallest `nameplate * dutycycle`
    then ratio = 1
    otherwise ratio = capacity / key.nameplate / key.dutycycle
    =#

    dict_effcap = Dict()
    for key in keys(selected_year_costs)
        effective_capacity = key.nameplate * key.dutycycle
        dict_effcap[effective_capacity] = key
    end
    list_effcap = sort(collect(keys(dict_effcap)))
    lowest_effcap = list_effcap[1]
    if (capacity < lowest_effcap)
       selected_key = dict_effcap[lowest_effcap]
       ratio = 1
    else
       effcap_less = findnearestleq(capacity, list_effcap)
       selected_key = dict_effcap[effcap_less]
       ratio = capacity / selected_key.nameplate / selected_key.dutycycle
    end
    selected_value = selected_year_costs[selected_key]

    result = selected_value.storage

    return(result)

end

function extrap_cost(x, d, fd)
   if x <= d[1]
      y = x * (fd[1] / d[1])
      if isnan(y)
         y = 0.0
      end
   else
      for i in 1:(length(d) - 1)
         if x > d[i] && x <= d[i + 1]
            y = fd[i] + (x - d[i]) * ((fd[i + 1] - fd[i]) / (d[i + 1] - d[i]))
         end
      end
   end

   return y
end

end
