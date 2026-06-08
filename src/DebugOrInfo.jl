#module T

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
macro debug_or_info(verbosity, main_message, keywords...)
    ekws = map(esc_rhs, keywords)
    quote
	if $(esc(verbosity)) > 0
            @info $(esc(main_message)) $(ekws...)
        else
            @debug $(esc(main_message)) $(ekws...)
        end
    end
end

function esc_rhs(eq)
    # assuming eq is Expr(:=, lhs, rhs)
    Expr(eq.head, eq.args[1], esc(eq.args[2]))
end


#end #module T
