"""
Miscellaneous functions.
"""
module Util

using Revise
# Export functions.
export parsebool
export parseenum

 string
"""
Base a boolean.
"""
function parsebool(text :: String)
    parseCandidates(Bool, [true, false], lowercase(text))
end


"""
Parse an enumeration.
"""
function parseenum(enum :: DataType, text)
    text = string(text)
    parseCandidates(enum, instances(enum), text)
end


"""
Parse a value from a list of candidates of a given type.
"""
function parseCandidates(dt, candidates, text)
    for candidate in candidates
        if text == string(candidate)
            return candidate
        end
    end
    error(string("Invalid ", dt, ": \"", text, "\"."))
end


end
