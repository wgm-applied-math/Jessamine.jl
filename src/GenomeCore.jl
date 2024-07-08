export AbstractGeneOp
export GenomeSpec, CellState, CellState, Instruction, Genome
export eval_time_step, op_eval, flat_workspace
export run_genome, num_instructions, num_operands, workspace_size
export short_show
export to_expr, CompiledGenome

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
    index_map::Vector{Tuple{Symbol,Int}}
end

function GenomeSpec(output_size, scratch_size, parameter_size, input_size, num_time_steps)
    index_map = vcat(
        [(:output, j) for j in 1:output_size],
        [(:scratch, j) for j in 1:scratch_size],
        [(:parameter, j) for j in 1:parameter_size],
        [(:input, j) for j in 1:input_size]
    )
    return GenomeSpec(output_size, scratch_size, parameter_size, input_size, num_time_steps, index_map)
end


"""
    workspace_size(g_spec::GenomeSpec)

Return the number of elements in the workspace vector specified by `g_spec`.
"""
function workspace_size(g_spec::GenomeSpec)
    return g_spec.output_size + g_spec.scratch_size + g_spec.parameter_size + g_spec.input_size
end

struct CellState{VOut,VScr,VPar,VIn}
    output::VOut
    scratch::VScr
    parameter::VPar
    input::VIn
    index_map::Vector{
        Tuple{
            Union{VOut,VScr,VPar,VIn},
            Int}}
end

function CellState(
    output::VOut,
    scratch::VScr,
    parameter::VPar,
    input::VIn)::CellState{VOut,VScr,VPar,VIn} where {VOut,VScr,VPar,VIn}
    V = Union{VOut,VScr,VPar,VIn}
    T = Tuple{V,Int}
    index_map = vcat(
        T[(output, j) for j in 1:length(output)],
        T[(scratch, j) for j in 1:length(scratch)],
        T[(parameter, j) for j in 1:length(parameter)],
        T[(input, j) for j in 1:length(input)]
    )
    return CellState{VOut,VScr,VPar,VIn}(output, scratch, parameter, input, index_map)
end


function Base.getindex(cs::CellState, i::Int)
    (v, j) = cs.index_map[i]
    return v[j]
end

function Base.getindex(cs::CellState, ix::AbstractVector)
    return [Base.getindex(cs, j) for j in ix]
end

function Base.length(cs::CellState)
    return length(cs.output) + length(cs.scratch) + length(cs.parameter) + length(cs.input)
end

function cell_output(cs::CellState{VOut,VScr,VPar,VIn})::VOut where {VOut,VScr,VPar,VIn}
    return cs.output
end

function flat_workspace(cs::CellState)
    return vcat(cs.output, cs.scratch, cs.parameter, cs.input)
end

@kwdef struct Instruction{OpType}
    op::OpType
    operand_ixs::Vector{Int}
end

function to_expr(g_spec::GenomeSpec, instr::Instruction, cs::Union{Symbol,Expr})
    field_fetches = [g_spec.index_map[k] for k in instr.operand_ixs]
    return to_expr(instr.op, cs, field_fetches)
end

function to_expr(g_spec::GenomeSpec, block::Vector{Instruction}, cs::Union{Symbol,Expr})
    if isempty(block)
        return :0
    elseif length(block) == 1
        return to_expr(g_spec, block[1], cs)
    else
        return Expr(:.+, (to_expr(g_spec, instr, cs) for instr in block)...)
    end
end

function to_expr(g_spec::GenomeSpec, blocks::Vector{Vector{Instruction}}, cs::Union{Symbol,Expr}, e_type::Union{Symbol,Expr})
    blocks_syn = [to_expr(g_spec, block, cs) for block in blocks]
    return :($e_type[ $(blocks_syn...) ])
end


"""
    op_eval(op::AbstractGeneOp, operands::AbstractVector)

Return `op` applied to a list of operands.
"""
function op_eval end

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
        function(cs::CellState{VOut,VScr,VPar,VIn}) where {VOut,VScr,VPar,VIn}
            new_output = $(to_expr(g_spec, genome.instruction_blocks[1:g_spec.output_size], :cs, :(eltype(VOut))))
            new_scratch = $(to_expr(g_spec, genome.instruction_blocks[scratch_begin:end], :cs, :(eltype(VScr))))
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
    num_operands(xs::AbstractVector)

Return the total number of instruction operands in `xs`.
"""
function num_operands(xs::AbstractVector)
    return sum(num_operands.(xs))
end


function eval_time_step(
    cell_state::CellState{Vector{Vector{T}},Vector{T},Vector{Vector{T}}},
    genome::Genome
    )::CellState{Vector{Vector{T}},Vector{T},Vector{Vector{T}}} where {T <: Number}

    num_rows = length(cell_state.input[1])

    new_output = map(genome.instruction_blocks[1:length(cell_state.output)]) do block
        val_out = zeros(T, num_rows)
        for instr in block
            val_out .= val_out .+ op_eval(instr.op, cell_state[instr.operand_ixs])
        end
        val_out
    end
    scratch_first = 1 + length(cell_state.output)
    new_scratch = map(genome.instruction_blocks[scratch_first:end]) do block
        val_scr = zeros(T, num_rows)
        for instr in block
            val_scr .= val_scr .+ op_eval(instr.op, cell_state[instr.operand_ixs])
        end
        val_scr
    end
    cell_state_next = CellState(
        new_output,
        new_scratch,
        cell_state.parameter,
        cell_state.input)
    return cell_state_next
end

function eval_time_step(
    cell_state::CellState{VOut,VScr,VPar,VIn},
    genome::Genome
    )::CellState{VOut,VScr,VPar,VIn} where {VOut,VScr,VPar,VIn}
    cell_state_next = CellState(
        deepcopy(cell_state.output),
        deepcopy(cell_state.scratch),
        cell_state.parameter,
        cell_state.input)
    for j in eachindex(genome.instruction_blocks)
        vec, i = cell_state_next.index_map[j]
        instructions = genome.instruction_blocks[j]
        val = zero_like(vec[i])
        for instr in instructions
            val = val .+ op_eval(instr.op, cell_state[instr.operand_ixs])
        end
        vec[i] = val
    end
    return cell_state_next
end


# scalar case
function zero_like(x::Number)
    return zero(x)
end

# vector case
function zero_like(v::AbstractVector)
    return zeros(eltype(v), size(v))
end

# scalar case
function zeros_like(x::Number, num_elts::Int)::Vector{typeof(x)}
    return zeros(typeof(x), num_elts)
end

# vector case
function zeros_like(v::T, num_elts::Int)::Vector{Vector{eltype(T)}} where {T <: AbstractVector}
    return [zeros(eltype(v), size(v)) for _ in 1:num_elts]
end

@kwdef struct CompiledGenome <: AbstractGenome
    genome::Genome
    expr::Expr
    f::Function

    function CompiledGenome(g_spec::GenomeSpec, genome::Genome)
        expr = to_expr(g_spec, genome)
        f = eval(expr)
        return new(genome, expr, f)
    end
end

function eval_time_step(
    cell_state::CellState{VOut,VScr,VPar,VIn},
    cg::CompiledGenome
    )::CellState{VOut,VScr,VPar,VIn} where {VOut,VScr,VPar,VIn}
    return invokelatest(cg.f, cell_state)
end

function num_operands(cg::CompiledGenome)
    return num_operands(cg.genome)
end

"""
    run_genome(g_spec::GenomeSpec, genome::Genome, parameter::AbstractVector, input::AbstractVector)::Vector

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
        parameter::VPar,
        input::VIn
    )::Vector{Vector{eltype(VIn)}} where {VPar,VIn}
    VOut = Vector{eltype(VIn)}
    @assert length(input) == g_spec.input_size
    @assert length(parameter) == g_spec.parameter_size
    output = zeros_like(input[1], g_spec.output_size)
    # output::VOut
    scratch = zeros_like(input[1], g_spec.scratch_size)
    # scratch::VOut
    current_state = CellState(output, scratch, parameter, input)
    outputs = Vector(undef, g_spec.num_time_steps)
    for t in 1:(g_spec.num_time_steps)
        future_state = eval_time_step(current_state, genome)
        # future_state::CellState{VOut,VOut,VPar,VIn}
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
