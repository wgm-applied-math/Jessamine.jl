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
macro debug_or_info(verbosity, messages...)
    quote
	if $(esc(verbosity)) > 0
            @info($(map(esc, messages)...))
        else
            @debug($(map(esc, messages)...))
        end
    end
end
