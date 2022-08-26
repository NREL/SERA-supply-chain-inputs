# Network Data Information

## Nodes
- Point locations where central production facilities can be located or where urban demand is known
- Demand locations and existing production locations are included
- nodes.tsv

## Links
- Links connect nodes (supply and/or demand locations)
- links.tsv

## Regions 
- Regions files must match with the price files so that correct prices can be assigned to production locations based on the node it is built at
- Balancing Areas - balancing-areas.tsv
- NERC (North American Electric Reliability Corporation) Regions - nerc-regions.tsv
  — With Canadian nodes: nerc-regions-nodes.tsv
	— Canadian nodes were assigned (wholly) to NERC regions: WECC-ExCA, MRO, and NPCC
	- NPCC: Ontario Quebec, Newfoundland, New Brunswick and Nova Scotia
	- MRO: Saskatchewan
	- WECC-ExCA: British Columbia and Alberta
- Census Divisions - census-divisions.tsv
- Territories - territories.tsv
- States - states.tsv
- USA - usa.tsv

## Existing Infrastructure
- Existing hydrogen production plants based on the 2018 IHS Chemical Handbook
- US and Canada locations are included
