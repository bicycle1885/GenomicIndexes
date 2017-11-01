# CSI
# ===
#
# Specification:
#   1: https://github.com/samtools/hts-specs/blob/master/CSIv1.pdf
#   2: https://github.com/samtools/hts-specs/blob/master/CSIv2.pdf
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

struct CSIBin
    # distinct bin
    bin::UInt32

    # virtual file offset
    loffset::UInt64

    # number of records stored in the bin
    n_rec::UInt64

    # chunks
    chunks::Vector{Chunk}
end

struct CSI
    # CSI version (0x01 or 0x02)
    version::UInt8

    # number of bits for the minimal interval
    minshift::Int32

    # depth of the binning index
    depth::Int32

    # auxilary data
    aux::Vector{UInt8}

    # bins
    bins::Vector{Vector{CSIBin}}

    # number of unmapped reads
    n_no_coor::Nullable{Int}
end

function Base.show(io::IO, csi::CSI)
    print(io, summary(csi), "(<version=$(csi.version),#seqs=$(length(csi.bins))>)")
end

"""
    CSI(input::AbstractString)

Load a CSI index from `input`.
"""
function CSI(input::AbstractString)
    return open(x->CSI(x), input)
end

"""
    CSI(input::IO)

Load a CSI index from `input`.
"""
function CSI(input_::IO)
    input = CodecBGZF.BGZFDecompressorStream(input_)

    # check magic bytes
    C = read(input, UInt8)
    S = read(input, UInt8)
    I = read(input, UInt8)
    v = read(input, UInt8)
    if C != UInt32('C') || S != UInt32('S') || I != UInt8('I') || v ∉ (0x01, 0x02)
        error("invalid csi magic bytes")
    end
    
    minshift = read(input, Int32)
    depth = read(input, Int32)
    l_aux = read(input, Int32)
    aux = read(input, l_aux)
    n_ref = read(input, Int32)
    bins = Vector{CSIBin}[]
    for _ in 1:n_ref
        n_bin = read(input, Int32)
        refbins = CSIBin[]
        for _ in 1:n_bin
            bin = read(input, UInt32)
            loffset = read(input, UInt64)
            if v == 0x02
                n_rec = read(input, UInt64)
            else
                n_rec = typemax(UInt64)
            end
            n_chunk = read(input, Int32)
            chunks = Chunk[]
            for _ in 1:n_chunk
                chunk_beg = read(input, UInt64)
                chunk_end = read(input, UInt64)
                push!(chunks, (chunk_beg, chunk_end))
            end
            push!(refbins, CSIBin(bin, loffset, n_rec, chunks))
        end
        push!(bins, refbins)
    end
    
    if !eof(input)
        n_no_coor = Nullable{Int}(read(input, UInt64))
    else
        n_no_coor = Nullable{Int}()
    end

    return CSI(v, minshift, depth, aux, bins, n_no_coor)
end

function overlapchunks(csi::CSI, seqid::Integer, interval::UnitRange{<:Integer})
    if !(1 ≤ seqid ≤ endof(csi.bins))
        throw(ArgumentError("sequence id $(seqid) is out of range"))
    elseif isempty(interval)
        return Chunk[]
    end
    chunks = Chunk[]
    bins = reg2bins(interval, csi.minshift, csi.depth)
    binindex = csi.bins[seqid]
    for bin in bins
        for x in binindex
            if x.bin == bin
                for chunk in x.chunks
                    if chunk[2] > x.loffset
                        push!(chunks, chunk)
                    end
                end
            end
        end
    end
    sort!(chunks)
    reduce!(chunks)
    return chunks
end

function reg2bins(interval, minshift, depth)
    return reg2bins(convert(UnitRange{Int}, interval), convert(Int, minshift), convert(Int, depth))
end

function reg2bins(interval::UnitRange{Int}, minshift::Int, depth::Int)
    bins = UInt32[]
    s = minshift + depth * 3
    bin_start = 0
    for d in 0:depth
        for k in ((first(interval)-1)>>s):((last(interval)-1)>>s)
            push!(bins, bin_start + k)
        end
        s -= 3
        bin_start = 8 * bin_start + 1
    end
    return bins
end
