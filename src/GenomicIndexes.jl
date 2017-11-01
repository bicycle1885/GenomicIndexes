module GenomicIndexes

import CodecBGZF

include("bgzfindex.jl")
include("bai.jl")
include("tabix.jl")
include("csi.jl")

end # module
