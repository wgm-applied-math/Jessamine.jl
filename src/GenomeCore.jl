export AbstractGeneOp
export GenomeSpec, CellState, Instruction, Genome
export eval_time_step, op_eval
export run_genome, num_instructions, num_operands, workspace_size
export short_show

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
end

"""
    workspace_size(g_spec::GenomeSpec)

Return the number of elements in the workspace vector specified by `g_spec`.
"""
function workspace_size(g_spec::GenomeSpec)
    return g_spec.output_size + g_spec.scratch_size + g_spec.parameter_size + g_spec.input_size
end

@kwdef struct CellState
    workspace::Vector
end

@kwdef struct Instruction{OpType}
    op::OpType
    operand_ixs::Vector{Int}
end

"""
    op_eval(op<:AbstractGeneOp, operands<:AbstractVector)

Return `op`` applied to a list of operands.
"""
function op_eval end

"""
    short_show([io::IO], x)

Print a short version of `x` to `io`, using `stdout` by default.
"""
short_show(x) = short_show(stdout, x)

""" A vector of blocks of instructions.  During a time step, each
instruction is evaluated on the current work space.  For each
`j`, the results of the instructions in block `j` are collected
and added, and this sum is used as the value of slot `j` in the
next work space.  """

@kwdef struct Genome
    instruction_blocks::Vector{Vector{Instruction}}
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
        cell_state::CellState,
        genome::Genome)
    workspace_next = deepcopy(cell_state.workspace)
    for dest in eachindex(genome.instruction_blocks)
        instructions = genome.instruction_blocks[dest]
        val = 0
        for instr in instructions
            val = val .+ op_eval(instr.op, cell_state.workspace[instr.operand_ixs])
        end
        workspace_next[dest] = val
    end
    return CellState(workspace_next)
end

"""
    run_genome(g_spec::GenomeSpec, genome::Genome, parameter::AbstractVector, input::AbstractVector)::Vector

Build a work space vector by concatenating zeros for each
instruction block in `genome`, followed by the `parameter`
vector, then the `input` vector.  Evaluate the instructions in
`genome`, repeating the evaluation `g_spec.num_time_steps`.
Return an array that contains, for each time step, the elements 1
through `g_spec.output_size` of the work space vector.
"""
function run_genome(
        g_spec::GenomeSpec,
        genome::Genome,
        parameter::AbstractVector,
        input::AbstractVector
)::Vector
    @assert length(input) == g_spec.input_size
    @assert length(parameter) == g_spec.parameter_size
    num_instr_blocks = length(genome.instruction_blocks)
    @assert num_instr_blocks == g_spec.output_size + g_spec.scratch_size
    scratch = zeros(num_instr_blocks)
    workspace = vcat(scratch, parameter, input)
    current_state = CellState(workspace)
    output_size = g_spec.output_size
    outputs = Vector(undef, g_spec.num_time_steps)
    for t in 1:(g_spec.num_time_steps)
        future_state = eval_time_step(current_state, genome)
        outputs[t] = future_state.workspace[1:output_size]
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
