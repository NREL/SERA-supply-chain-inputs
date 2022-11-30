# FIXME: This file should probably be named something else.


import SERA.Supply.Types: MaterialName, MaterialUsage, Prices, ZoneID


"""
Look up and interpolate prices.

FIXME: document arguments and add example.
"""
function pricelookup(prices :: Prices, usage :: MaterialUsage, material :: MaterialName, zone :: ZoneID, year :: Int16, quantity :: Float64, interpolate=false) :: Float64
    # First,  find how much material is already used in the zone and year.
    # Second, the new quantity used will be that plus the amount being priced.
    # Third, find the two points on the price curve that bracket the year.
    # Fourth, either interpolate or just take the first value.
end


"""
Return a function for looking up prices.

FIXME: document arguments and add example.
"""
function pricefunction(prices :: Prices, usage :: MaterialUsage, material :: MaterialName, zone :: ZoneID, interpolate=false)
    function f(year :: Int16, quantity :: Float64)
        pricelookup(prices, usage, material, zone, interpolate)
    end
    return f
end
