export recombine

"""
    recombine(rng, g_left, g_right)

Given two genomes, line up the instruction blocks.  For each
pair, choose a prefix of the block from the left, drop its length
from the block from the right, and concatenate the resulting
arrays to form a new block.

The resulting recombined blocks are used to build a new `Genome`.
No mutations are applied.
"""

function recombine(rng::AbstractRNG, g_left::Genome, g_right::Genome)
    @assert length(g_left.instruction_blocks) == length(g_right.instruction_blocks)
    mixed = similar(g_left.instruction_blocks)
    for j in eachindex(g_left.instruction_blocks)
        l_block = g_left.instruction_blocks[j]
        r_block = g_right.instruction_blocks[j]
        l_m = length(l_block)
        if l_m == 0
            l_k = 1
            l_frag = []
        elseif l_m == 1
            l_k = 1
            l_frag = l_block
        else
            l_k = rand(rng, DiscreteUniform(1, l_m - 1))
            l_frag = l_block[1:l_k]
        end
        r_frag = r_block[(l_k + 1):end]
        new_block = vcat(l_frag, r_frag)
        mixed[j] = new_block
    end
    return Genome(mixed)
end
