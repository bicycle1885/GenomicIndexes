module GenomicIndexes

export
    BAI,
    Tabix,
    CSI,
    seqnames

import CodecBGZF

include("bgzfindex.jl")
include("bai.jl")
include("tabix.jl")
include("csi.jl")

end # module
