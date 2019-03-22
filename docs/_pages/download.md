---
layout: page
title: Download
permalink: /download/
order: 4
share: false
---

- TOC
{:toc}

## Index

{% for item in site.data.download-index %}
### {{ item.organism }}
  {% for data in item.data %}
<li>{{ data[0] }}</li>
<table style="border-collapse: collapse; border: none;">
{% for genome in data[1] %}
<tr style="border: none;"><td style="border: none;">{{ genome[0] }}</td><td style="border: none;"><a href="{{ genome[1].url }}">{{ genome[1].url }}</a></td></tr>
{% endfor %}
</table>
{% endfor %}
{% endfor %}


    genome: HISAT2 index for reference
    genome_snp: HISAT2 Graph index for reference plus SNPs
    genome_tran: HISAT2 Graph index for reference plus transcripts
    genome_snp_tran: HISAT2 Graph index for reference plus SNPs and transcripts

## Binaries

