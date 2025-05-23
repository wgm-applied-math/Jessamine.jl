# Make language server happy
if false
end

export AbstractGeneOp
export GenomeSpec, CellState, Instruction
export AbstractGenome, Genome
export v_convert, v_unconvert
export eval_time_step, op_eval, flat_workspace
export run_genome, num_instructions, num_operands, workspace_size
export short_show
export to_expr, compile, CompiledGenome
export zeros_like, zero_like

"""
An operation to be performed as part of the function of a gene.
"""
abstract type AbstractGeneOp end

"""
A collection of parameters specifying the genome architecture and mutation processes.
"""
@kwdef struct GenomeSpec
    "How many output slots are in the state array"
    output_size::Int

    "How many scratch slots are in the state array"
    scratch_size::Int

    "How many parameter slots are in the state array"
    parameter_size::Int

    "How many input slots are in the state array"
    input_size::Int

    "How many time steps in an evaluation cycle"
    num_time_steps::Int

    "Table of which field and index therein corresponds to an overall index"
    index_map::Vector{Tuple{Symbol, Int}}
end

function GenomeSpec(output_size, scratch_size, parameter_size, input_size, num_time_steps)
    index_map = vcat(
        [(:output, j) for j in 1:output_size],
        [(:scratch, j) for j in 1:scratch_size],
        [(:parameter, j) for j in 1:parameter_size],
        [(:input, j) for j in 1:input_size]
    )
    return GenomeSpec(
        output_size, scratch_size, parameter_size, input_size, num_time_steps, index_map)
end

function Base.convert(::Type{GenomeSpec}, g_spec::GenomeSpec)
    return g_spec
end

function Base.convert(::Type{GenomeSpec}, x)
    return GenomeSpec(
        x.output_size,
        x.scratch_size,
        x.parameter_size,
        x.input_size,
        x.num_time_steps)
end

"""
    workspace_size(g_spec::GenomeSpec)

Return the number of elements in the workspace vector specified by `g_spec`.
"""
function workspace_size(g_spec::GenomeSpec)
    return g_spec.output_size + g_spec.scratch_size + g_spec.parameter_size +
           g_spec.input_size
end

struct CellState{E} <: AbstractVector{E}
    output::Vector{E}
    scratch::Vector{E}
    parameter::Vector{E}
    input::Vector{E}
    index_map::Vector{Tuple{Vector{E}, Int}}
end

function v_convert(::Type{<:AbstractArray}, x::AbstractArray)
    return x
end

function v_convert(ref::Type{<:AbstractArray}, x)
    singleton = [convert(eltype(ref), x)]
    return convert(ref, singleton)
end

function v_convert(ref::Type, x)
    return convert(ref, x)
end

function v_unconvert(xv::AbstractArray)
    @assert length(xv) == 1
    return xv[1]
end

function CellState(
        output,
        scratch,
        parameter,
        input)
    # E = least common supertype of the element types of all of these vectors
    E = typejoin(eltype(output), eltype(scratch), eltype(parameter), eltype(input))
    n_output = length(output)
    n_scratch = length(scratch)
    n_parameter = length(parameter)
    n_input = length(input)
    n_total = n_output + n_scratch + n_parameter + n_input
    index_map = Vector{Tuple{Vector{E}, Int}}(undef, n_total)
    for j in 1:n_output
        index_map[j] = (output, j)
    end
    offset = n_output
    for j in 1:n_scratch
        index_map[offset + j] = (scratch, j)
    end
    offset += n_scratch
    for j in 1:n_parameter
        index_map[offset + j] = (parameter, j)
    end
    offset += n_parameter
    for j in 1:n_input
        index_map[offset + j] = (input, j)
    end
    return CellState{E}(output, scratch, parameter, input, index_map)
end

function Base.getindex(cs::CellState, i::Int)
    (v, j) = cs.index_map[i]
    return v[j]
end

function Base.getindex(cs::CellState, ix::AbstractArray)
    return [Base.getindex(cs, j) for j in ix]
end

function Base.length(cs::CellState)
    return length(cs.output) + length(cs.scratch) + length(cs.parameter) + length(cs.input)
end

function Base.size(cs::CellState)
    return [length(cs)]
end

function cell_output(cs::CellState)
    return cs.output
end

function flat_workspace(cs::CellState)
    return vcat(cs.output, cs.scratch, cs.parameter, cs.input)
end

@kwdef struct Instruction
    op::AbstractGeneOp
    operand_ixs::Vector{Int}
end

# *sigh*
function Instruction(op, operands::Vector{Any})
    @assert isempty(operands)

    return Instruction(op, Int[])
end

function to_expr(g_spec::GenomeSpec, instr::Instruction, cs::Union{Symbol, Expr})
    field_fetches = [g_spec.index_map[k] for k in instr.operand_ixs]
    return to_expr(instr.op, cs, field_fetches)
end

function to_expr(g_spec::GenomeSpec, block::Vector{Instruction}, cs::Union{Symbol, Expr})
    if isempty(block)
        return :[0.0]
    elseif length(block) == 1
        return to_expr(g_spec, block[1], cs)
    else
        return Expr(:call, :.+, (to_expr(g_spec, instr, cs) for instr in block)...)
    end
end

function to_expr(g_spec::GenomeSpec, blocks::Vector{Vector{Instruction}},
        cs::Union{Symbol, Expr}, e_type::Union{Symbol, Expr})
    blocks_syn = [to_expr(g_spec, block, cs) for block in blocks]
    return :($e_type[$(blocks_syn...)])
end

"""
    op_eval(op::AbstractGeneOp, operands::AbstractArray)

Return `op` applied to a list of operands.
"""
function op_eval end

"""
    op_eval_add_into!(dest, op, workspace, indices)

Apply the operator `op` to `operands = workspace[indices]` and add the result elementwise
in place to dest.  The default is to do
`dest .= dest .+ op_eval(op, workspace, indices)`.

The idea here is that in-place accumulation instead of the
default implementation can sometimes avoid the allocation of a
temporary array to hold the result of `op_eval(...)`.

"""
function op_eval_add_into!(
        dest::AbstractArray, op::AbstractGeneOp, workspace::AbstractArray, indices::AbstractArray)
    dest .+= op_eval(op, workspace, indices)
end

"""
    short_show([io::IO], x)

Print a short version of `x` to `io`, using `stdout` by default.
"""
short_show(x) = short_show(stdout, x)

"""Abstract base type for genomes."""
abstract type AbstractGenome end

""" A vector of blocks of instructions.  During a time step, each
instruction is evaluated on the current work space.  For each
`j`, the results of the instructions in block `j` are collected
and added, and this sum is used as the value of slot `j` in the
next work space.  """
@kwdef struct Genome <: AbstractGenome
    instruction_blocks::Vector{Vector{Instruction}}
end

function to_expr(g_spec::GenomeSpec, genome::Genome)
    scratch_begin = 1 + g_spec.output_size
    quote
        function (cs::CellState{E}) where {E}
            new_output = $(to_expr(
                g_spec, genome.instruction_blocks[1:(g_spec.output_size)], :cs, :E))
            new_scratch = $(to_expr(
                g_spec, genome.instruction_blocks[scratch_begin:end], :cs, :E))
            return CellState(new_output, new_scratch, cs.parameter, cs.input)
        end
    end
end

"""
    num_instructions(genome::Genome)

Return the total number of instructions in all blocks in the genome.
"""
function num_instructions(genome::Genome)
    return sum(length.(genome.instruction_blocks))
end

"""
    num_operands(genome::Genome)

Return the total number of operands in all instructions in the genome.
"""
function num_operands(genome::Genome)
    return num_operands(genome.instruction_blocks)
end

"""
    num_operands(instruction::Instruction)

Return the total number of operands in the instruction.
"""
function num_operands(instruction::Instruction)
    return length(instruction.operand_ixs)
end

"""
    num_operands(xs::AbstractArray)

Return the total number of instruction operands in `xs`.
"""
function num_operands(xs::AbstractArray)
    return sum(num_operands.(xs))
end

function eval_time_step(
        cell_state::CellState{E},
        genome::Genome
)::CellState{E} where {E <: AbstractArray}
    # local new_output::Vector{E}
    n_output = length(cell_state.output)
    new_output = [zero_like(cell_state.input[1]) for j in 1:n_output]
    for j in 1:n_output
        block = genome.instruction_blocks[j]
        for instr in block
            op_eval_add_into!(new_output[j], instr.op, cell_state, instr.operand_ixs)
        end
    end
    n_scratch = length(cell_state.scratch)
    new_scratch = [zero_like(cell_state.input[1]) for j in 1:n_scratch]
    for j in 1:n_scratch
        block = genome.instruction_blocks[n_output + j]
        for instr in block
            op_eval_add_into!(new_scratch[j], instr.op, cell_state, instr.operand_ixs)
        end
    end
    cell_state_next = CellState(
        new_output,
        new_scratch,
        cell_state.parameter,
        cell_state.input)
    return cell_state_next
end

function eval_time_step(
        cell_state::CellState{E},
        genome::Genome
)::CellState{E} where {E <: Number}
    # local new_output::Vector{E}
    new_output = zero_like(cell_state.output)
    n_output = length(new_output)
    for j in 1:n_output
        block = genome.instruction_blocks[j]
        for instr in block
            new_output[j] += op_eval(instr.op, cell_state, instr.operand_ixs)
        end
    end
    scratch_first = 1 + length(cell_state.output)
    # local new_scratch::Vector{E}
    new_scratch = zero_like(cell_state.scratch)
    n_scratch = length(new_scratch)
    for j in 1:length(new_scratch)
        block = genome.instruction_blocks[n_output + j]
        for instr in block
            new_scratch[j] += op_eval(instr.op, cell_state, instr.operand_ixs)
        end
    end
    cell_state_next = CellState(
        new_output,
        new_scratch,
        cell_state.parameter,
        cell_state.input)
    return cell_state_next
end

function zero_like(ref::T)::T where {T <: Number}
    return convert(T, 0)
end

# vector case
# This is a fallback implementation that works for GPU arrays and other oddities
function zero_like(v::V)::AbstractArray where {V <: AbstractArray}
    return v .* 0
end

function zero_like(v::Array)::Array
    return zeros(eltype(v), size(v))
end

function zeros_like(
        v::AbstractArray, num_elts::Int)
    return [zero_like(v) for _ in 1:num_elts]
end

function zeros_like(v::T, num_elts::Int)::Vector{T} where {T <: Number}
    return zeros(typeof(v), num_elts)
end

@kwdef struct CompiledGenome <: AbstractGenome
    genome::Genome
    expr::Expr
    f::Function
end

function compile(g_spec::GenomeSpec, genome::Genome)
    expr = to_expr(g_spec, genome)
    try
        f = eval(expr)
        return CompiledGenome(genome, expr, f)
    catch
        short_show(genome)
        @show expr
        rethrow()
    end
end

function compile(g_spec::GenomeSpec, cg::CompiledGenome)
    return cg
end

function eval_time_step(
        cell_state::CellState{E},
        cg::CompiledGenome
)::CellState{E} where {E}
    return invokelatest(cg.f, cell_state)
end

function num_operands(cg::CompiledGenome)
    return num_operands(cg.genome)
end

"""
    run_genome(g_spec::GenomeSpec, genome::Genome, parameter::AbstractArray, input)::Vector

Build a work space vector using zeros for each output and scratch
slot, followed by the `parameter` vector, then the `input`
vector.  Evaluate the instructions in `genome`, repeating the
evaluation `g_spec.num_time_steps`.  Return an array that
contains, for each time step, the elements 1 through
`g_spec.output_size` of the work space vector.

"""
function run_genome(
        g_spec::GenomeSpec,
        genome::AbstractGenome,
        parameter::AbstractArray,
        input
)
    input_v = separate_columns(input)
    @assert length(input_v) == g_spec.input_size
    @assert length(parameter) == g_spec.parameter_size
    parameter_v = v_convert.(eltype(input_v), parameter)
    output_v = zeros_like(input_v[1], g_spec.output_size)
    scratch_v = zeros_like(input_v[1], g_spec.scratch_size)
    current_state = CellState(output_v, scratch_v, parameter_v, input_v)
    outputs = Vector(undef, g_spec.num_time_steps)
    for t in 1:(g_spec.num_time_steps)
        future_state = eval_time_step(current_state, genome)
        outputs[t] = cell_output(future_state)
        current_state = future_state
    end
    return outputs
end

function short_show(io::IO, g::Genome)
    for dest in eachindex(g.instruction_blocks)
        block = g.instruction_blocks[dest]
        print(io, "$(dest) = sum [ ")
        for instr in block
            print(io, "")
            short_show(io, instr.op)
            for j in instr.operand_ixs
                print(io, " $j")
            end
            print(io, ", ")
        end
        println(io, "]")
    end
end

function short_show(io::IO, cg::CompiledGenome)
    short_show(io, cg.genome)
    println(io, cg.expr)
end

function check_genome(g)
    for block in g.instruction_blocks
        for instr in block
            if isempty(instr.operand_ixs)
                return false
            end
        end
    end
    return true
end
