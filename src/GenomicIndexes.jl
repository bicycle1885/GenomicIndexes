module GenomicIndexes

export
    BAI,
    Tabix,
    CSI

import CodecBGZF

include("bgzfindex.jl")
include("bai.jl")
include("tabix.jl")
include("csi.jl")

end # module
