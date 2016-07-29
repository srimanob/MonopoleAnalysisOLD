# MonopoleAnalysis
How to checkout the package for analysis
--------------
Please use the fresh CMSSW, not recycle from the ntuple maker.
<pre><code>cmsrel CMSSW_8_0_13
cd CMSSW_8_0_13/src
cmsenv
git init
git clone -b 80X https://github.com/srimanob/MonopoleAnalysis
cd MonopoleAnalysis/NtupleAnalyzer
scram b
</code></pre>

If you want to checkout specific tag
git checkout tags/80X_pre5 [No need in general. Just checkout the head one in 80X]

How to run analysis code
--------------


How to commit new version of code
--------------
After you checkout the head, and make modification
<pre><code>git checkout 80X
git add [ files that you edited ]
git commit -m 'comment'
git push origin 80X
git tag -a [version] -m 'comment' [check existing tags first]
git push --tags
</code></pre>
