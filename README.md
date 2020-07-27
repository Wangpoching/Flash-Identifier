# Flash_Identifier

Flash_Identifier is a toolkit for extracting flash signals from a timestamp data.
#### Liscence: MIT <http://choosealicense.com/licenses/mit/>

## Installation & Initialization

### Requirements 

* Linux
* R 3.4 and higher

`$ chmod -rwx Flash_Identifier.R`
## Test
Download a timestamp data measured from *Abscondita  cerata*:
```
$ wget https://github.com/Kinggerm/GetOrganelleGallery/raw/master/Test/reads/Arabidopsis_simulated.1.fq.gz
```
Then do the fast flash identify:

`$ ./Flash_Identifier.R -180411NK_BM1.csv`

You are going to get the same result as here.
## Instruction
#### Preparing Data
Currently, Flash_Identifier.R was written for CSV file. Please note that the first column recored the timestamp and the second column recorded the measured value. 


#### Flash Identifier Flowchart
![GITHUB](https://drive.google.com/file/d/13VaI8umcaArPHzBy7pGhuUs0-xYfYFR5/view?usp=sharing "Flash Identifier Flowchart")


## Recipes
To identify *Luciola kagiana* flash, I use:

`$ ./Flah_Identifier.R -180413NK_RF5.csv`

or I want to open the debug mode:

`$ ./Flah_Identifier.R -180413NK_RF5.csv --debug`

If I want to define an environmental background 10 (defalt=0) by myself:

`$ ./Flah_Identifier.R -180413NK_RF5.csv --debug --en 10`

or just use the median value of the valleys as environmental background:

`$ ./Flah_Identifier.R -180413NK_RF5.csv --debug --en_median`

or just use the minimum value of the valleys as environmental background:

`$ ./Flah_Identifier.R -180413NK_RF5.csv --debug --en_minimum`

If I want to set a peak filter coefficient 1 (defalt=0.3) by myself:

`$ ./Flah_Identifier.R -180413NK_RF5.csv --debug --filter 1`



## Contact
If your question is running specific, please attach the flash_ide.log.txt file
* Report bugs & Open issues here.
* Send email to me ([b04208021@g.ntu.edu.tw](https://mail.google.com/mail/u/0/?view=cm&fs=1&tf=1&source=mailto&to=b04208021@g.ntu.edu.tw))








