---
layout: compress
# WARNING: Don't use '//' to comment out code, use '{% comment %}' and '{% endcomment %}' instead.
---

{%- comment -%}
  See: <https://docs.mathjax.org/en/latest/options/input/tex.html#tex-options>
{%- endcomment -%}

MathJax = {
  tex: {
    {%- comment -%} start/end delimiter pairs for in-line math {%- endcomment -%}
    inlineMath: [
      ['$', '$'],
      ['\\(', '\\)']
    ],
    {%- comment -%} start/end delimiter pairs for display math {%- endcomment -%}
    displayMath: [
      ['$$', '$$'],
      ['\\[', '\\]']
    ],
    {%- comment -%} equation numbering: 'none' or 'ams' or 'all' {%- endcomment -%}
    tags: 'all',
    {%- comment -%} custom macros used across the posts {%- endcomment -%}
    macros: {
      Re:  ['\\operatorname{Re}\\left({#1}\\right)', 1],
      Im:  ['\\operatorname{Im}\\left({#1}\\right)', 1],
      EV:  ['\\mathbb{E}\\left({#1}\\right)', 1],
      bs:  ['\\boldsymbol{#1}', 1],
      vc:  ['\\mathrm{vec}\\left({#1}\\right)', 1],
      cov: ['\\mathrm{cov}\\left({#1}\\right)', 1],
      var: ['\\mathrm{var}\\left({#1}\\right)', 1],
      std: ['\\mathrm{std}\\left({#1}\\right)', 1],
      overbar: ['\\mkern 1.5mu\\overline{\\mkern-1.5mu#1\\mkern-1.5mu}\\mkern 1.5mu', 1]
    }
  }
};
