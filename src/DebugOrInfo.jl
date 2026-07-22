# module T

export @debug_or_info

"""
    @debug_or_info(verbosity, messages...)

If the expression `verbosity` evaluates to a number greater than
0, call `@info` with `messsges`, otherwise call `@debug` with
messages.

That is, the messges are used to create a `@debug` log item,
unless verbosity is turned on, in which case it's promoted to a
`@info` log item.
"""
macro debug_or_info(verbosity, exs...)
    # This is adapted directly from julia/base/logging/logging.jl
    Base.CoreLogging.logmsg_code(
        (Base.CoreLogging.@_sourceinfo)...,
        :($(esc(verbosity)) > 0 ? Base.CoreLogging.Info : Base.CoreLogging.Debug),
        exs...)
end

# end #module T
