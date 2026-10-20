#!/usr/bin/env bash

prefetch SRR3734796
prefetch SRR3734797
prefetch SRR3734798
prefetch SRR3734820
prefetch SRR3734824
prefetch SRR3734825

fasterq-dump SRR3734796 -O SRR3734796
fasterq-dump SRR3734797 -O SRR3734797   
fasterq-dump SRR3734798 -O SRR3734798
fasterq-dump SRR3734820 -O SRR3734820
fasterq-dump SRR3734824 -O SRR3734824
fasterq-dump SRR3734825 -O SRR3734825