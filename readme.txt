Overview
--------
Droog is a python script for generating misspelling lists for trade and generic drug
names.  It automates a four-step pipeline process:
- generate all 1 and 2 edit distance variants of the correct spelling
- compute misspelling probability based on a noisy channel model
- pass the highest probability variants through a Google domain filter
- select the number of high page count variants required by the application

Author
------
Robert D. Hogan, Ph.D.
Terminologix LLC
700 Fifth Avenue, Suite 3A
Antigo, WI 54424
rhogan@terminologix.com

License
-------
Droog is made available under a GPL v3 license.  Contact author for other licensing 
options including closed commercial use.

Files
-----
/droog
	/data
		misspellings.txt			# list of misspellings for noisy channel model
		char_counts.txt				# character counts used in noisy channel model
		bigram_counts.txt			# bigram counts used in noisy channel model
	/results
		/drug1						# result folder for each drug analyzed
			drug1_d1_variants.txt	# edit distance 1 variants of drug name with probabilities
			drug1_d2_variants.txt	# edit distance 2 variants of drug name with probabilities
			drug1_likely_variants.txt	# highest probability variants
			drug1_filtered_variants.txt	# high probability variants with positive Google hits
		/drug2
	__init__.py
	droog.py						# CLI interface
	noisy.py						# Droog noisy channel model calculations
	configuration.py				# user configurations

The *_d*_variants.txt files contain comma delimited lists with the following fields
- variant text
- edit distance (1 or 2)
- probability of first edit
- probability of second edit
- edit(s) (e.g. ted-mrd = deletion of 'e' after 't' then deletion of 'r' after 'm')

Requirements
------------
- Python 3.5
- matplotlib
- Google Custom Search Engine (CSE) business account

Setup
-----
1. download and install files in folder "droog"
2. edit configuration.py to add CSE search engine ID and key
3. setup search engine to limit search to relevant domains
3. (optionally) modify configuration constants to select input to Google filter

Usage
-----
runs from command line
usage: droog.py [-h] -d DRUG [-g] [-f] [-a] [-l LIST] [-p]

Drug misspelling generator

optional arguments:
  -h, --help            show this help message and exit
  -d DRUG, --drug DRUG  correctly spelled drug name
  -g, --generate        generate variants with probabilities
  -f, --filter          filter variants with Google CSE
  -a, --analyze         print table of terms vs. cumulative counts
  -l LIST, --list LIST  list n, n% or 'all' filtered terms
  -p, --plot            plot cumulative page count vs. list number
  
run options in order g > f > [a | l | p] for a given drug name

Notes
-----
1.	Google query string must have double quotes around it to prevent pre-processing
2.	noisy.py contains all logic for noisy channel model for use in non-CLI program
3.	noisy channel model uses add-one smoothing with lower limit for bigram counts
	set to 250 to avoid spurious high probability calculations.
4.	Google search engine may be pointed to any set of relevant domains.  We recommend
	dailystrength.org, inspire.com and healingwell.com/community for analyzing drug
	names.  Do not use the entire web for the search.
5.	Google results may include false positive variants that are unrelated to the
	drug.  For example, the variant "acteria" of the drug "actemra" is actually a deletion
	edit of "bacteria".  It is necessary to scan the result list to check for ambiguous
	variants.  Variants can be tested using the Google web interface with a query of the form:
	["variant" site:dailystrength.org OR site:inspire.com OR site:healingwell.com/community]
	Square brackets are not included.  Variant name must be in quotes to avoid query
	pre-processing by Google search engine.  Placing a "#" in front of the name in the 
	drug_filtered_variants.txt file causes that variant to be ignored by the "analyze",
	"list", and "plot" command line options.