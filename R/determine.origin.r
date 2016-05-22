
    
  determine.origin = function(id) {
    
    # determine study id
    
    study=NULL
    if (id %in% c(0:600))      study=1
    if (id %in% c(2350:2399))  study=5
    if (id %in% c(1000:1600))  study=6
    if (id %in% c(2401:2403, 2411:2450)) study=7
    if (id %in% c(1601:2349))  study=8
    if (id %in% c(6000:6349))  study=9
    if (id %in% c(2716:2850))  study=10
    if (id %in% c(2456:2715))  study=11
    if (id %in% c(5050:5099))  study=12
    if (id %in% c(5100:5149))  study=13
    if (id %in% c(2851:2900, 3051:3262, 3446:3545))  study=14
    if (id %in% c(3263:3444, 3546:3600))  study=15
    if (id %in% c(4000:4246))  study=16
    if (id %in% c(7450:7542))  study=17
    if (id %in% c(7543:7699))  study=18
    if (id %in% c(4798:4999, 5250:5520, 6350:6999))  study=19
      
    if (id %in% c(5521:5808))  study=20
    if (id %in% c(4250:4480))  study=21
    if (id %in% c(4482:4797))  study=22
    if (id %in% c(7000:7449, 7700:7999))  study=23
    if (id %in% c(5809:5999, 8501:8570, 9229:9649))  study=24
    if (id %in% c(8571:8828))  study=25
    if (id %in% c(8829:9228))  study=26
    if (id %in% c(10482:10499, 11082:11334, 11350:11387))  study=27
    if (id %in% c(4301, 8040:8050, 8137:8150, 8201:8298, 9684:9741, 9748:9784, 10254:10406, 10468:10481, 11039:11081 ))  study=28
    if (id %in% c(4248:4249, 8000:8039, 8051:8068, 8070:8136, 8151:8200, 10407:10467, 11000:11036))  study=29
    if (id %in% c(8299:8500, 9650:9683))  study=30
    if (id %in% c(10032:10253))  study=31
    if (id %in% c(9742:9747, 9785:10031))  study=32
    if (id %in% c(10501:10530, 10609:10716))  study=33
    if (id %in% c(10531:10608, 10718:10798))  study=34
    if (id %in% c(10800:10949, 11400:11449))  study=35
    if (id %in% c(10950:10999, 11335:11349, 11388:11399, 11450:11499, 13000:13149))  study=36
    if (id %in% c(13150:13490))  study=37
    if (id %in% c(13491:13649, 15650:15673))  study=38
    if (id %in% c(13650:13834))  study=39
    if (id %in% c(13835:14109))  study=40
    if (id %in% c(14110:14226))  study=41
    if (id %in% c(14227:14306, 14308:14311, 14313:14314, 14316, 14318:14345, 14347:14430))  study=42
    if (id %in% c(14431:14880))  study=43
    if (id %in% c(15750:15779, 15800:15899, 15950:15999))  study=44
    if (id %in% c(15780:15799, 15900:15949))  study=45

    if (id %in% c(14881:14951))  study=46
    if (id %in% c(14952:15001))  study=47
    if (id %in% c(15002:15199))  study=48
    if (id %in% c(15400:15499, paste("t", 1605:1659, sep="")))  study=49
    if (id %in% c(15500:15649, 15674:15749, 16000:16711))  study=50
    if (id %in% c(16712:17199, 17300:17399, 17500:17607))  study=51

    if (id %in% c(15200:15399, paste("t", c(1217:1350, 1601:1603), sep="") ))  study=52
    if (id %in% c(17608:17849, 17953))  study=53
    if (id %in% c(17200:17299, 17850:17999))  study=54
    if (id %in% c(18043:18099))  study=55
    if (id %in% c(17401:17499, 18000:18042, 18100:18158,  paste("t", 1676:1694, sep="")))  study=56
    if (id %in% c(18226:18326, paste("t", 3026:3175, sep=""), paste("s", 99102:99196, sep="")))  study=57

      return (study)
    
  }
   


