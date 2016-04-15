import threading

when isMainModule:
  threadsInit()
  echo threadNum, "/", numThreads
  threads:
    echo threadNum, "/", numThreads
