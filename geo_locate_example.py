<html>
<head>
<title>geo_locate_example.py</title>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
<style type="text/css">
.s0 { color: #cf8e6d;}
.s1 { color: #bcbec4;}
.s2 { color: #bcbec4;}
.s3 { color: #7a7e85;}
.s4 { color: #6aab73;}
.s5 { color: #2aacb8;}
</style>
</head>
<body bgcolor="#1e1f22">
<table CELLSPACING=0 CELLPADDING=5 COLS=1 WIDTH="100%" BGCOLOR="#606060" >
<tr><td><center>
<font face="Arial, Helvetica" color="#000000">
geo_locate_example.py</font>
</center></td></tr></table>
<pre><span class="s0">import </span><span class="s1">numpy </span><span class="s0">as </span><span class="s1">np</span>
<span class="s0">import </span><span class="s1">matplotlib</span><span class="s2">.</span><span class="s1">pyplot </span><span class="s0">as </span><span class="s1">plt</span>
<span class="s0">from </span><span class="s1">scipy </span><span class="s0">import </span><span class="s1">constants</span><span class="s2">, </span><span class="s1">optimize</span>
<span class="s0">from </span><span class="s1">pyproj </span><span class="s0">import </span><span class="s1">Geod</span>
<span class="s0">import </span><span class="s1">pickle</span>
<span class="s0">import </span><span class="s1">pandas </span><span class="s0">as </span><span class="s1">pd</span>
<span class="s0">from </span><span class="s1">scipy</span><span class="s2">.</span><span class="s1">optimize </span><span class="s0">import </span><span class="s1">LinearConstraint</span>
<span class="s0">from </span><span class="s1">scipy</span><span class="s2">.</span><span class="s1">optimize </span><span class="s0">import </span><span class="s1">minimize</span>
<span class="s3"># Load the node config dictionary this contains information about the nodes</span>
<span class="s3"># (receivers) position and name</span>
<span class="s1">node_config </span><span class="s2">= </span><span class="s1">np</span><span class="s2">.</span><span class="s1">load</span><span class="s2">(</span><span class="s4">&quot;node_config_fixed.npy&quot;</span><span class="s2">, </span><span class="s1">allow_pickle</span><span class="s2">=</span><span class="s0">True</span><span class="s2">).</span><span class="s1">item</span><span class="s2">()</span>

<span class="s3"># Now load in waveforms</span>
<span class="s1">waveforms </span><span class="s2">= </span><span class="s1">np</span><span class="s2">.</span><span class="s1">load</span><span class="s2">(</span><span class="s4">&quot;exported_waveforms.npy&quot;</span><span class="s2">, </span><span class="s1">allow_pickle</span><span class="s2">=</span><span class="s0">True</span><span class="s2">)</span>

<span class="s3"># For each waveform the following variables are available</span>
<span class="s3"># envelope      : The magnitude of the waveform envelope with time</span>
<span class="s3"># FT            : The Fast fourier transform of samples</span>
<span class="s3"># FT_recon      : The Fast fourire transform of the reconstructed time series</span>
<span class="s3"># node_unique_id: A 17 character unique identifier for a node</span>
<span class="s3"># peak_index    : The index in the raw data stream that the peak in envelope resides</span>
<span class="s3"># prod_time     : The time the waveform was extracted</span>
<span class="s3"># raw_data      : The raw time series the waveform was extracted from that has been high pass filtered at 1kHz</span>
<span class="s3"># recon_data    : The reconstructed complex time series the waveform was extracted from same as reconstructed samples</span>
<span class="s3"># samples       : The real part of extracted waveform</span>
<span class="s3"># site_name     : Name of site waveform was recorded at</span>
<span class="s3"># t1            : The sum of samples ^ 2</span>
<span class="s3"># tdv_penalty   : A time difference variance penalty </span>
<span class="s3"># timestamp     : np.datetime64 object which contains the timestamp at index 256 of samples, envelope, raw_data and recon_data</span>

<span class="s3"># The sampling resolution of LEELA is 9.142 us</span>



<span class="s1">node_code </span><span class="s2">= {}</span>
<span class="s0">for </span><span class="s1">key </span><span class="s0">in </span><span class="s1">node_config</span><span class="s2">.</span><span class="s1">keys</span><span class="s2">():</span>
    <span class="s1">node_code</span><span class="s2">[</span><span class="s1">node_config</span><span class="s2">[</span><span class="s1">key</span><span class="s2">][</span><span class="s4">'Site'</span><span class="s2">]] = </span><span class="s1">key</span>



<span class="s4">''' 
This function finds the Theoretical Arrival Time Difference. This is done through calculating 
the distances between the guess and each node. distGA is the distance between the guess and  
the reference node. distGBs is a list of all the distances between the guess and all other nodes. 
This is done using the geod function in pyproj which calculates distances between 2 points taking 
into account the geodisic lines of the earth (WGS84 being the best representation of the earth). 
Theoretical Arrival Time Difference is calculated using both distGA and distGB as well as the phase 
velocity (vp) which will be a little over the speed of light. The list TATD is outputted to be utilised 
in the RES function. 
'''</span>

<span class="s0">def </span><span class="s1">TATD</span><span class="s2">(</span><span class="s1">waveforms</span><span class="s2">, </span><span class="s1">guess</span><span class="s2">):</span>
    <span class="s1">dist_GBs </span><span class="s2">= []</span>
    <span class="s1">geod </span><span class="s2">= </span><span class="s1">Geod</span><span class="s2">(</span><span class="s1">ellps</span><span class="s2">=</span><span class="s4">'WGS84'</span><span class="s2">)</span>
    <span class="s1">_</span><span class="s2">, </span><span class="s1">_</span><span class="s2">, </span><span class="s1">dist_GA </span><span class="s2">= </span><span class="s1">geod</span><span class="s2">.</span><span class="s1">inv</span><span class="s2">(</span><span class="s1">waveforms</span><span class="s2">[</span><span class="s5">0</span><span class="s2">][</span><span class="s4">'lon'</span><span class="s2">], </span><span class="s1">waveforms</span><span class="s2">[</span><span class="s5">0</span><span class="s2">][</span><span class="s4">'lat'</span><span class="s2">], </span><span class="s1">guess</span><span class="s2">[</span><span class="s5">0</span><span class="s2">], </span><span class="s1">guess</span><span class="s2">[</span><span class="s5">1</span><span class="s2">])</span>

    <span class="s0">for </span><span class="s1">i </span><span class="s0">in </span><span class="s1">range</span><span class="s2">(</span><span class="s5">1</span><span class="s2">, </span><span class="s5">7</span><span class="s2">):</span>
        <span class="s1">_</span><span class="s2">, </span><span class="s1">_</span><span class="s2">, </span><span class="s1">dist_GB </span><span class="s2">= </span><span class="s1">geod</span><span class="s2">.</span><span class="s1">inv</span><span class="s2">(</span><span class="s1">waveforms</span><span class="s2">[</span><span class="s1">i</span><span class="s2">][</span><span class="s4">'lon'</span><span class="s2">], </span><span class="s1">waveforms</span><span class="s2">[</span><span class="s1">i</span><span class="s2">][</span><span class="s4">'lat'</span><span class="s2">], </span><span class="s1">guess</span><span class="s2">[</span><span class="s5">0</span><span class="s2">], </span><span class="s1">guess</span><span class="s2">[</span><span class="s5">1</span><span class="s2">])</span>
        <span class="s1">dist_GBs</span><span class="s2">.</span><span class="s1">append</span><span class="s2">(</span><span class="s1">dist_GB</span><span class="s2">)</span>
        <span class="s3"># print(waveforms[x]['site_name'], waveforms[x]['lat'], waveforms[x]['lon'], dist_GB)</span>
    <span class="s1">TATDs </span><span class="s2">= []</span>

    <span class="s1">vp </span><span class="s2">= </span><span class="s5">2.99792458e8 </span><span class="s2">* </span><span class="s5">1.00425</span>

    <span class="s0">for </span><span class="s1">i </span><span class="s0">in </span><span class="s1">range</span><span class="s2">(</span><span class="s5">0</span><span class="s2">, </span><span class="s5">6</span><span class="s2">):</span>
        <span class="s1">TATD </span><span class="s2">= </span><span class="s1">abs</span><span class="s2">(</span><span class="s1">dist_GBs</span><span class="s2">[</span><span class="s1">i</span><span class="s2">] / </span><span class="s1">vp</span><span class="s2">) - </span><span class="s1">abs</span><span class="s2">(</span><span class="s1">dist_GA </span><span class="s2">/ </span><span class="s1">vp</span><span class="s2">)</span>
        <span class="s1">TATDs</span><span class="s2">.</span><span class="s1">append</span><span class="s2">(</span><span class="s1">TATD</span><span class="s2">)</span>
    <span class="s0">return </span><span class="s1">TATDs</span>



<span class="s4">''' 
RES function takes in an intial guess of where the lightening point will be, 
the waveforms dictionary containing data necessary to find the TATD (see TATD 
function) and OATDs is the list of OATDs (see above). The guess will change after 
each iteration of the minimisation function is complete. The point of this 
fuction is to calculate the residual (res) of how close the guess is to the actual 
lightening strike location. The smaller res gets, the closer the guess is to the 
location of the lightening strike. 
'''</span>

<span class="s0">def </span><span class="s1">RES</span><span class="s2">(</span><span class="s1">guess</span><span class="s2">, </span><span class="s1">waveforms</span><span class="s2">, </span><span class="s1">OATDs</span><span class="s2">):</span>
    <span class="s3"># Finds TATDs</span>
    <span class="s1">TATDs </span><span class="s2">= </span><span class="s1">TATD</span><span class="s2">(</span><span class="s1">waveforms</span><span class="s2">, </span><span class="s1">guess</span><span class="s2">)</span>
    <span class="s3"># Finds variance</span>
    <span class="s1">vari </span><span class="s2">= </span><span class="s1">np</span><span class="s2">.</span><span class="s1">var</span><span class="s2">(</span><span class="s1">OATDs</span><span class="s2">)</span>
    <span class="s3"># determines N from amount of values in OATDs</span>
    <span class="s1">N </span><span class="s2">= </span><span class="s1">len</span><span class="s2">(</span><span class="s1">OATDs</span><span class="s2">)</span>

    <span class="s3"># Calculate the value of RES</span>
    <span class="s1">sumATD </span><span class="s2">= </span><span class="s5">0</span>
    <span class="s0">for </span><span class="s1">i </span><span class="s0">in </span><span class="s1">range</span><span class="s2">(</span><span class="s5">0</span><span class="s2">, </span><span class="s5">5</span><span class="s2">):</span>
        <span class="s1">sumATD </span><span class="s2">= </span><span class="s1">sumATD </span><span class="s2">+ (</span><span class="s1">TATDs</span><span class="s2">[</span><span class="s1">i</span><span class="s2">] - </span><span class="s1">OATDs</span><span class="s2">[</span><span class="s1">i</span><span class="s2">]) ** </span><span class="s5">2 </span><span class="s2">/ </span><span class="s1">vari</span>
    <span class="s1">RES </span><span class="s2">= </span><span class="s1">np</span><span class="s2">.</span><span class="s1">sqrt</span><span class="s2">((</span><span class="s5">1 </span><span class="s2">/ (</span><span class="s1">N </span><span class="s2">- </span><span class="s5">2</span><span class="s2">) * </span><span class="s1">sumATD</span><span class="s2">))</span>
    <span class="s1">print</span><span class="s2">(</span><span class="s1">RES</span><span class="s2">)</span>
    <span class="s0">return </span><span class="s1">RES</span>


<span class="s4">''' 
This adds longitude and and latitude to the waveforms dictionary for ease to only need to pass  
through one dictionary to all functions. A list is also defined to hold all observed arrival 
time differences. This is the time difference between our reference node (waveforms[0]) and  
another arbitrary node. This is obviously doesn't vary so it makes sense to not attach this 
to a function. This reference node is consistent when calculating distances between  
nodes as well as the theoretical arrival time difference. 
'''</span>
<span class="s1">OATDs </span><span class="s2">= []</span>
<span class="s0">for </span><span class="s1">i </span><span class="s0">in </span><span class="s1">range</span><span class="s2">(</span><span class="s5">0</span><span class="s2">,</span><span class="s5">7</span><span class="s2">):</span>
    <span class="s1">waveforms</span><span class="s2">[</span><span class="s1">i</span><span class="s2">][</span><span class="s4">'lat'</span><span class="s2">] = </span><span class="s1">node_config</span><span class="s2">[</span><span class="s1">node_code</span><span class="s2">[</span><span class="s1">waveforms</span><span class="s2">[</span><span class="s1">i</span><span class="s2">][</span><span class="s4">'site_name'</span><span class="s2">]]][</span><span class="s4">'Position'</span><span class="s2">][</span><span class="s4">'lat'</span><span class="s2">]</span>
    <span class="s1">waveforms</span><span class="s2">[</span><span class="s1">i</span><span class="s2">][</span><span class="s4">'lon'</span><span class="s2">] = </span><span class="s1">node_config</span><span class="s2">[</span><span class="s1">node_code</span><span class="s2">[</span><span class="s1">waveforms</span><span class="s2">[</span><span class="s1">i</span><span class="s2">][</span><span class="s4">'site_name'</span><span class="s2">]]][</span><span class="s4">'Position'</span><span class="s2">][</span><span class="s4">'lon'</span><span class="s2">]</span>
<span class="s0">for </span><span class="s1">i </span><span class="s0">in </span><span class="s1">range</span><span class="s2">(</span><span class="s5">1</span><span class="s2">, </span><span class="s5">7</span><span class="s2">):</span>
    <span class="s1">OATDs</span><span class="s2">.</span><span class="s1">append</span><span class="s2">((</span><span class="s1">waveforms</span><span class="s2">[</span><span class="s1">i</span><span class="s2">][</span><span class="s4">'timestamp'</span><span class="s2">] - </span><span class="s1">waveforms</span><span class="s2">[</span><span class="s5">0</span><span class="s2">][</span><span class="s4">'timestamp'</span><span class="s2">]).</span><span class="s1">astype</span><span class="s2">(</span><span class="s1">float</span><span class="s2">)/</span><span class="s5">1e9</span><span class="s2">)</span>
<span class="s0">for </span><span class="s1">i </span><span class="s0">in </span><span class="s1">range</span><span class="s2">(</span><span class="s5">0</span><span class="s2">, </span><span class="s5">6</span><span class="s2">):</span>
    <span class="s1">print</span><span class="s2">(</span><span class="s1">OATDs</span><span class="s2">[</span><span class="s1">i</span><span class="s2">])</span>


<span class="s3">#Initalise guess in london</span>
<span class="s4">''' 
guess = [0, 0] 
guess[0] = -0.118092 #longitude 
guess[1] = 51.509865 #latitude 
'''</span>
<span class="s1">guess </span><span class="s2">= </span><span class="s1">np</span><span class="s2">.</span><span class="s1">array</span><span class="s2">([-</span><span class="s5">0.118092</span><span class="s2">, </span><span class="s5">51.509865</span><span class="s2">])</span>
<span class="s4">''' 
This calls the minimise function in scipy which will minimise RES function until RES is sufficenitly  
small. Minimise allows the use of a variety of diiferent methods to minimise but we are using the Nelder-Mead 
method. 
'''</span>
<span class="s1">result </span><span class="s2">= </span><span class="s1">minimize</span><span class="s2">(</span><span class="s1">RES</span><span class="s2">, </span><span class="s1">guess</span><span class="s2">, </span><span class="s1">args</span><span class="s2">=(</span><span class="s1">waveforms</span><span class="s2">, </span><span class="s1">OATDs</span><span class="s2">), </span><span class="s1">method</span><span class="s2">=</span><span class="s4">'Nelder-Mead'</span><span class="s2">)</span>



<span class="s1">print</span><span class="s2">(</span><span class="s1">result</span><span class="s2">)</span>

</pre>
</body>
</html>