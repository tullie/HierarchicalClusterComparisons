Cluster comparison
==================

How to run?
-----------
Within the build.gradle file:
  - Change the 'mainClassName' field to the algorithm driver of choice (BirchDriver, ChameleonDriver, ChameleonExtensionDriver).
  - Change the 'args' variable. The first value expects the input file name - see example files for format. The second args value contains the number of clusters to be found.

Run:
``gradle run``

Resulting labeled files will be saved in src/main/resources/. These can be viewed with:
``python vis/plot.py <filepath>``


TODO
----
- Output for CURE and ROCK (don't try to run these yet - working on it now).
- Work out how to specify cluster size for BIRCH.
- Finish Chameleon extension implementation.
