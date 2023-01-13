# Cache-Simulator

Author: Sailasya Gangavarapu

Designed and developed a flexible cache and memory hierarchy simulator in Python

CACHE is configurable in terms of supporting any cache size, associativity, and block size, specified at the beginning of simulation:

  1. SIZE: Total bytes of data storage.
  2. ASSOC: The associativity of the cache (ASSOC = 1 is a direct-mapped cache).
  3. BLOCKSIZE: The number of bytes in a block.

Successfully implemented three different Replacement Policies namely:

  1. LRU (Least Recently Used)
  2. FIFO (First In First Out)
  3. Optimal Policy

Inputs to the Simulator:

  1. BLOCKSIZE
  2. L1_SIZE
  3. L1_ASSOC
  4. L2_SIZE
  5. L2_ASSOC
  6. REPLACEMENT_POLICY (0 for LRU, 1 for FIFO, 2 for Optimal)
  7. INCLUSION_PROPERTY (0 for non-inclusive, 1 for inclusive)
  8. trace_file

Outputs from the Simulator:

  1. number of L1 reads
  2. number of L1 read misses
  3. number of L1 writes
  4. number of L1 write misses
  5. L1 miss rate
  6. number of writebacks from L1 to next level
  7. number of L2 reads
  8. number of L2 read misses if there is a L2 cache
  9. number of L2 writes
  10. number of L2 write misses
  11. L2 miss rate
  12. number of writebacks from L2 to memory
  13. total memory traffic

Commands to run the code:

- Open command prompt and navigate to the folder with the sim_cache.py code.

- Run the following command once in the code directory:

    python sim_cache.py BLOCKSIZE L1_SIZE L1_ASSOC L2_SIZE L2_ASSOC REPLACEMENT_POLICY INCLUSION_PROPERTY trace_file
