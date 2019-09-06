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

{% assign target = site.data.download-binary.latest_version %}
{% for release in site.data.download-binary.release %}
{% assign version = release['version'] %}
{% if version == target or target == null %}
{% assign name = release['name'] %}
### Version: {{name}} {{version}}
<table style="border-collapse: collapse; border: none;">
<tr style="border: none;"><td style="border: none;" colspan="2"><b>Release Date</b>: {{release['date']}}</td></tr>
{% for artifact in release['artifacts'] %}
{% assign type = artifact[0] %}
<tr style="border: none;"><td style="border: none;">{{type}}</td><td style="border: none;"><a href="{{artifact[1]}}">{{artifact[1]}}</a></td></tr>
{% endfor %}
</table>
{% endif %}
{% endfor %}

