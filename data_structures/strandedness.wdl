version 1.3

## Possible strandedness protocols used during RNA-Seq library prep.
enum Strandedness {
    ## The protocol is unknown or otherwise unspecified.
    ##
    ## Some tooling can automatically derive a suspected strandedness protocol,
    ## if one is not specified.
    Unspecified,
    ## The RNA-Seq library was prepped with an unstranded protocol.
    Unstranded,
    ## The RNA-Seq library was prepped with a Stranded-Reverse protocol.
    StrandedReverse,
    ## The RNA-Seq library was prepped with a Stranded-Forward protocol.
    StrandedForward,
}
