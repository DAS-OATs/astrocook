from astrocook.cookbook_absorbers import CookbookAbsorbers as Absorbers
from astrocook.cookbook_continuum import CookbookContinuum as Continuum
from astrocook.cookbook_edit import CookbookEdit as Edit
from astrocook.cookbook_flux import CookbookFlux as Flux
from astrocook.cookbook_general import CookbookGeneral as General
from astrocook.cookbook_synthetic import CookbookSynthetic as Synthetic
from astrocook.cookbook_templates import CookbookTemplates as Templates
from astrocook.cookbook_view import CookbookView as View
import inspect
import sys

classes = [General, Flux, Continuum, Absorbers, Synthetic, Templates]
names = [c.__name__ for c in classes]
if sys.argv[1]=='Other':
    ind = 6
    cs = [Edit, View]
else:
    ind = names.index('Cookbook'+sys.argv[1])
    cs = [classes[ind]]
ms = sys.argv[2].split(',')

preamble = \
'---\n\
layout: default\n\
title: %s cookbook\n\
parent: Cookbooks\n\
nav_order: %i\n\
---\n\
\n\
# %s cookbook\n\
{: .no_toc}\n\
\n\
%s\n\
\n\
\n\
## Table of contents\n\
{: .no_toc .text-delta }\n\
\n\
1. TOC\n\
{:toc}\n\
---\n\
'

table = \
'<table>\n\
  <tbody>\n\
    <tr>\n\
      <td style="vertical-align:top;width:200px"><strong>Method</strong></td>\n\
      <td style="vertical-align:top"><code>%s</code></td>\n\
    </tr>\n\
    <tr>\n\
      <td style="vertical-align:top"><strong>Parameters</strong></td>\n\
      <td style="vertical-align:top">\n\
%s\
      </td>\n\
    </tr>\n\
    <tr>\n\
      <td style="vertical-align:top;width:200px"><strong>JSON template</strong></td>\n\
      <td style="vertical-align:top"><pre>\n\
%s\
    </pre></td>\n\
    </tr>\n\
  </tbody>\n\
</table>'

for c in cs:
    i = c()
    if c == cs[0]:
        doc = str(i.__doc__)
        details = doc.split("@details")[-1]
        details = ''.join(details.split('\n'))
        details = ' '.join(details.split())
        print(preamble % (sys.argv[1],ind,sys.argv[1],details))
    for m in ms:
        if '#' in m:
            print(' '.join(m.split('_')))
        elif hasattr(i, m):
            f = getattr(i, m)
            doc = str(f.__doc__)

            brief = doc\
                        .split("@brief")[-1]\
                        .split("@details")[0]
            print("### %s" % brief)

            pars = doc\
                      .split("@return")[0]\
                      .split("@param")[1:]
            pars = [' '.join(p.split()) for p in pars]
            #print(doc)
            #print(pars)

            mstr = c.__name__+'.'+m
            parstr = ''
            if len(pars)==0:
                parstr += '        –\n'
            else:
                parstr += '        <ul>\n'
                for p in pars:
                    k = p.split(' ')[0]
                    v = ' '.join(p.split(' ')[1:])
                    parstr += '          '
                    parstr += '<li><code>%s</code>: %s</li>\n' % (k,v)
                parstr += '        </ul>\n'
            sig = str(inspect.signature(f))
            jstr = \
                '{\n'\
                '  "cookbook": "cb",\n'\
                '  "recipe": "%s",\n'\
                '  "params": {\n' % m
            if len(pars) > 0:
                for s in sig[1:-1].split(', '):
                    k = s.split('=')[0]
                    v = s.split('=')[-1].split('(')[-1].split(')')[0]
                    v = ''.join(v.split('"'))
                    if v == "True": v = "true"
                    if v == "False": v = "false"
                    if v == "None": v = "null"
                    jstr += '    "%s": "%s",\n' % (k,v)
                jstr = jstr[:-2]+'\n'
            jstr += \
                '  }\n'\
                '}'
            print(table % (mstr, parstr, jstr))
            print()

            details = doc\
                         .split("@details")[-1]\
                         .split("@param")[0]\
                         .split("@return")[0]
            details = ''.join(details.split('\n'))
            details = ' '.join(details.split())
            details1 = details.split("@")
            details2 = details1[-1].split(' ')
            details = '`'.join([details1[0], details2[0], ' '+' '.join(details2[1:])])
            print("%s" % details)
        print()
