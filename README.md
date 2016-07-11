# MonopoleAnalysis
How to checkout the package for analysis
--------------
Please use the fresh CMSSW, not recycle from the ntuple maker.
<pre><code>cmsrel CMSSW_8_0_13
cd CMSSW_8_0_13
cmsenv
git init
git clone -b 80X https://github.com/srimanob/MonopoleAnalysis
cd MonopoleAnalysis
git checkout tags/80X_pre2
scram b
cd MonojetAnalysis/NtupleAnalyzer
</code></pre>

How to run analysis code
--------------


How to commit new version of code
--------------
<pre><code>After "git checkout"
git checkout 80X
git add [ files that you edited ]
git commit -m 'comment'
git push origin 80X
git tag -a [version] -m 'comment'
git push --tags
</code></pre>
