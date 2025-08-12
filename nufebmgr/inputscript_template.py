#jinja2 template for inputscript
#using as a python string instead of as a j2 file for ease of package managment

TEMPLATE_STR = \
"""{%- for line in header -%}
{{ line }}
{% endfor %}

{% for section in system_settings -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{% if values | length < 2 %}{{ "\n" }}# {{ values[0] }} {% else %}{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}{% endif %}
{% endfor %}
{%- endfor %}

{%- for section in microbes_and_groups %}
{{ section.title }}
    {%- for entry in section.bug_groups %}
{{ entry.values() | list | join(' ')}}
    {%- endfor %}
    {%- for entry in section.neighbors %}
{{ entry.values() | list | join(' ')}}
    {%- endfor %}
{% endfor %}

{%- for section in mesh_grid_and_substrates %}
{{ section.title }}
    {%- for entry in section.content %}
{{ entry.values() | list | join(' ')}}
    {%- endfor %}
{% endfor %}

{%- for section in biological_processes %}
  {%- for subsection, items in section.items() %}
    {%- if subsection == 'title' %}
{{ items }}
    {%- else %}
      {%- for entry in items %}
        {%- set values = entry.values() | list %}
        {%- if values | length < 2 %}

# {{ values[0] }}
        {%- else %}
{{ values | join(' ') }}
        {%- endif %}
      {%- endfor %}
    {%- endif %}
  {%- endfor %}
{%- endfor %}



{% for section in physical_processes -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{%- if values | length < 2 -%}
{{ "\n" }}# {{ values[0] }} 
{%- else -%}
{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}
{% endif %}
{% endfor %}
{% endfor %}

{%- for section in post_physical_processes %}
{{ section.title }}
    {%- for item in section.boundary_layer %}
        {%- if 'title' in item %}
#{{ item.title }}
         {%- else %}
{{ item.name }} {{item.fix_name}} {{item.group}} {{item.fix_loc}} {{item.distance_m}} 1 zhi {% if item.comment %} # {{ item.comment }}{% endif %}
        {%- endif %}
    {%- endfor %}
{% endfor %}

{%- for section in chemical_processes %}
{{ section.title }}

    {%- for item in section.diffusion_biofilm_ratios %}
        {%- if 'title' in item %}

#{{ item.title }}
         {%- else %}
{{ item.name }} {{item.fix_name}} {{item.group}} {{item.fix_loc}} {{item.sub1}} ratio {{item.coeff1}} {% if item.comment %} # {{ item.comment }}{% endif %}
        {%- endif %}
    {%- endfor %}
    
    {%- for item in section.diffusion_coefficients %}
        {%- if 'title' in item %}
#{{ item.title }}
         {%- else %}
{{ item.name }} {{item.fix_name}} {{item.group}} {{item.fix_loc}} {{item.sub1}} {{item.coeff1}} {% if item.comment %} # {{ item.comment }}{% endif %}
        {%- endif %}
    {%- endfor %}
{% endfor %}

{%- for section in computation_output %}
  {%- for subsection, items in section.items() %}
    {%- if subsection == 'title' %}
{{ items }}
    {%- else %}
      {%- for entry in items %}
        {%- set values = entry.values() | list %}
        {%- if values | length < 2 %}

# {{ values[0] }}
        {%- else %}
{{ values | join(' ') }}
        {%- endif %}
      {%- endfor %}
    {%- endif %}
  {%- endfor %}
{%- endfor %}

{% for section in run -%}
{{ section.title }}
{% for entry in section.content -%}
{%- set values = entry.values() | list -%}
{%- if values | length < 2 -%}
{{ "\n" }}# {{ values[0] }} 
{%- else -%}
{{ values[0] }}\t\t{{ values[1:-1] | join(' ') }}{% if values[-1] != '' %}\t\t#{{ values[-1] }}{% endif -%}
{% endif %}
{% endfor %}
{% endfor %}

"""