#!/usr/bin
#
# @author James Boocok

sed "s|!BASE_PIPELINE!|`pwd`|g"  config/default_template.yaml > pipelines/defaults.yaml 

