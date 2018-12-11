``` r
library(microbiomefhir)

library(tidyverse)
library(curatedMetagenomicData)
library(phyloseq)
library(here)
```

The microbiomfhir package allows the conversion of a patient’s OTU table
into a structured FHIR DiagnosticReport. At this point the package only
supports the import of [phyloseq](https://joey711.github.io/phyloseq/)
objects.

Planned features
================

-   integrate SNOMED API calls to code for bacterial species
-   integrate PROCEDURE ressource to annotate more experiment related
    information

Installation Instructions
=========================

The package is currenlty available via GitHub

``` r
devtools::install_github("https://github.com/NiklasTR/microbiomefhir")
```

The package has, among others, the following requirements:

-   tidyverse toolbox
-   phyloseq
-   curatedMetagenomicData

Example
=======

Let’s use the package’s built-in function to collect an example from
Levi Waldron’s excellent curated Metagenomic Data
[ressource](http://waldronlab.io/curatedMetagenomicData/).

``` r
pseq <- load_example()
```

Here is a glimpse at one of the OTU table profiles. We can see the
patient identifier as the column name. Each row represents a taxonomic
unit (from kingdom - Bacteria down to many species).

``` r
otu_table(pseq)[,1] %>% head()
```

    ## OTU Table:          [6 taxa and 1 samples]
    ##                      taxa are rows
    ##                   LomanNJ_2013.metaphlan_bugs_list.stool:OBK1122
    ## k__Bacteria                                            100.00000
    ## p__Bacteroidetes                                        83.55662
    ## p__Firmicutes                                           15.14819
    ## p__Proteobacteria                                        1.29520
    ## c__Bacteroidia                                          83.55662
    ## c__Clostridia                                           14.93141

We can simply export the patient’s OTU table into a FHIR
DiagnosticReport:

``` r
parse_json(pseq, file = here("output/my_first_microbiome_fhir_report.json"))
```

We can take a closer look at out output JSON file.

``` r
readLines(here("output/my_first_microbiome_fhir_report.json"), 100)
```

    ##   [1] "{"                                                                         
    ##   [2] "  \"resourceType\": [\"DiagnosticReport\"],"                               
    ##   [3] "  \"id\": [\"101\"],"                                                      
    ##   [4] "  \"meta\": {"                                                             
    ##   [5] "    \"tag\": ["                                                            
    ##   [6] "      {"                                                                   
    ##   [7] "        \"system\": \"http://example.org/fhir/CodeSystem/workflow-codes\","
    ##   [8] "        \"code\": \"01\","                                                 
    ##   [9] "        \"display\": \"Needs Review\""                                     
    ##  [10] "      }"                                                                   
    ##  [11] "    ]"                                                                     
    ##  [12] "  },"                                                                      
    ##  [13] "  \"text\": {"                                                             
    ##  [14] "    \"status\": [\"generated\"],"                                          
    ##  [15] "    \"div\": [\"<div Lorem Ipsum; ????<\\/p><\\/div>\"]"                   
    ##  [16] "  },"                                                                      
    ##  [17] "  \"contained\": {"                                                        
    ##  [18] "    \"resourceType\": [\"Specimen\"],"                                     
    ##  [19] "    \"id\": [\"stool\"],"                                                  
    ##  [20] "    \"accessionIdentifier\": ["                                            
    ##  [21] "      {"                                                                   
    ##  [22] "        \"system\": \"http://example.org/specimens/2019\","                
    ##  [23] "        \"value\": \"X352356\""                                            
    ##  [24] "      }"                                                                   
    ##  [25] "    ],"                                                                    
    ##  [26] "    \"status\": [\"unavailable\"],"                                        
    ##  [27] "    \"type\": ["                                                           
    ##  [28] "      {"                                                                   
    ##  [29] "        \"coding\": ["                                                     
    ##  [30] "          {"                                                               
    ##  [31] "            \"system\": \"http://snomed.info/sct\","                       
    ##  [32] "            \"code\": \"119339001\","                                      
    ##  [33] "            \"display\": \"Stool specimen\""                               
    ##  [34] "          }"                                                               
    ##  [35] "        ]"                                                                 
    ##  [36] "      }"                                                                   
    ##  [37] "    ],"                                                                    
    ##  [38] "    \"subject\": ["                                                        
    ##  [39] "      {"                                                                   
    ##  [40] "        \"reference\": \"Patient/example\""                                
    ##  [41] "      }"                                                                   
    ##  [42] "    ],"                                                                    
    ##  [43] "    \"receivedTime\": [\"2019-08-16T07:03:00Z\"],"                         
    ##  [44] "    \"collection\": ["                                                     
    ##  [45] "      {"                                                                   
    ##  [46] "        \"collector\": {"                                                  
    ##  [47] "          \"display\": \"Patient\""                                        
    ##  [48] "        },"                                                                
    ##  [49] "        \"collectedDateTime\": \"2011-08-16T06:15:00Z\","                  
    ##  [50] "        \"method\": {"                                                     
    ##  [51] "          \"coding\": ["                                                   
    ##  [52] "            {"                                                             
    ##  [53] "              \"system\": \"http://hl7.org/fhir/v2/0488\","                
    ##  [54] "              \"code\": \"SWA\""                                           
    ##  [55] "            }"                                                             
    ##  [56] "          ]"                                                               
    ##  [57] "        }"                                                                 
    ##  [58] "      }"                                                                   
    ##  [59] "    ],"                                                                    
    ##  [60] "    \"issued\": [null],"                                                   
    ##  [61] "    \"valueCodeableConcept\": ["                                           
    ##  [62] "      {"                                                                   
    ##  [63] "        \"coding\": {}"                                                    
    ##  [64] "      }"                                                                   
    ##  [65] "    ],"                                                                    
    ##  [66] "    \"valueQuantity\": ["                                                  
    ##  [67] "      {}"                                                                  
    ##  [68] "    ],"                                                                    
    ##  [69] "    \"12\": ["                                                             
    ##  [70] "      {"                                                                   
    ##  [71] "        \"resourceType\": \"Observation\","                                
    ##  [72] "        \"id\": \"k__Bacteria\","                                          
    ##  [73] "        \"status\": \"final\","                                            
    ##  [74] "        \"subject\": {"                                                    
    ##  [75] "          \"display\": \"Roel\""                                           
    ##  [76] "        },"                                                                
    ##  [77] "        \"issued\": \"2019-08-16 07:03:00\","                              
    ##  [78] "        \"valueCodeableConcept\": {"                                       
    ##  [79] "          \"coding\": ["                                                   
    ##  [80] "            {"                                                             
    ##  [81] "              \"system\": \"http://snomed.info/sct\","                     
    ##  [82] "              \"display\": \"Bacteria\""                                   
    ##  [83] "            }"                                                             
    ##  [84] "          ]"                                                               
    ##  [85] "        },"                                                                
    ##  [86] "        \"valueQuantity\": {"                                              
    ##  [87] "          \"value\": 100,"                                                 
    ##  [88] "          \"unit\": \"percent\","                                          
    ##  [89] "          \"code\": \"%\""                                                 
    ##  [90] "        }"                                                                 
    ##  [91] "      }"                                                                   
    ##  [92] "    ],"                                                                    
    ##  [93] "    \"13\": ["                                                             
    ##  [94] "      {"                                                                   
    ##  [95] "        \"resourceType\": \"Observation\","                                
    ##  [96] "        \"id\": \"p__Bacteroidetes\","                                     
    ##  [97] "        \"status\": \"final\","                                            
    ##  [98] "        \"subject\": {"                                                    
    ##  [99] "          \"display\": \"Roel\""                                           
    ## [100] "        },"
