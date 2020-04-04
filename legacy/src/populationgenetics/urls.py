from django.conf.urls import patterns, include, url
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf import settings
# Uncomment the next two lines to enable the admin:
# from django.contrib import admin
# admin.autodiscover()

urlpatterns = patterns('simulations.views',
    # Examples:
    # url(r'^$', 'definettiweb.views.home', name='home'),
    # url(r'^definettiweb/', include('definettiweb.foo.urls')),

    # Uncomment the admin/doc line below to enable admin documentation:
    # url(r'^admin/doc/', include('django.contrib.admindocs.urls')),

    # Uncomment the next line to enable the admin:
    # url(r'^admin/', include(admin.site.urls)),
	url(r'^$','index'),
    url(r'^simulations/$','index'),
    url(r'^simulations/new_index/$','new_index'),
    url(r'^simulations/definetti/$','definetti'),
    url(r'^simulations/definetti/plot.png$','plot'),
    url(r'^simulations/definetti/run/$','run'),
    url(r'^simulations/mutation/$','mutation'),
    url(r'^simulations/mutation/plot.png$','plot_mutation'),
    url(r'^simulations/mutation/run/$','run_mutation'),
    url(r'^simulations/finiteh/$','finiteh'),
    url(r'^simulations/finiteh/run/$', 'run_finiteh'),
    url(r'^simulations/finiteh/plot.png$','plot_finiteh'),
    url(r'^simulations/mullers/$','mullers'),
    url(r'^simulations/mullers/run/$', 'run_mullers'),
    url(r'^simulations/mullers/plot.png$','plot_mullers'),
    url(r'^simulations/naturalselection/$','naturalselection'),
    url(r'^simulations/naturalselection/run/$', 'run_naturalselection'),
    url(r'^simulations/naturalselection/plot.png$','plot_naturalselection'),
    url(r'^simulations/qtl/$','qtl'),
    url(r'^simulations/qtl/run/$', 'run_qtl'),
    url(r'^simulations/qtl/plot.png$','plot_qtl'),
    url(r'^simulations/driftselection/$','driftselection'),
    url(r'^simulations/driftselection/run/$', 'run_driftselection'),
    url(r'^simulations/driftselection/plot.png$','plot_driftselection'),
	url(r'^simulations/neutraltheory/$','neutral_theory'),
	url(r'^simulations/neutraltheory/run/$', 'run_neutral_theory'),
	url(r'^simulations/neutraltheory/plot.png$','plot_neutral_theory'),
    url(r'^simulations/nucsubstitution/$','nuc_substitution'),
    url(r'^simulations/nucsubstitution/run/$', 'run_nuc_substitution'),
    url(r'^simulations/nucsubstitution/plot.png$','plot_nuc_substitution'),

)

#if settings.DEBUG:
#    urlpatterns += patterns('django.contrib.staticfiles.views',
#                            url(r'^static/(?P<path>.*)$','serve'),
#                            {'document_root': '/Users/marieidowu/pythonCode/definettiweb/definettiPlotter/static/'},
#                            )