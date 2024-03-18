library(mutscan)

output <- digestFastqs(
  fastqForward = "/home/ubuntu/AWTallyNN/tallynn/python/240318_IDTGblock_Flipped_2",
  elementsForward = "SPVS",
  primerForward = "GCTGGTGAGGTTGCGGATAACG",
  elementLengthsForward = c(-1, 22, 18, -1),
  maxReadLength = 2000  
  # Add other arguments as needed
)