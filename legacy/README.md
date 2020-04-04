population-genetics
===================

Setting up application on test server for cython version
Using pyximport in views.py

Need to create .pyxbld directory under /var/www owned by apache

https://github.com/cython/cython/blob/master/pyximport/pyximport.py

VirtualHost
------------
mod_wsgi
The application uses third party C Extension modules and this can cause the app to crash.
"
Another major cause of crashes is third party C extension modules for Python which are not compatible with being used in a Python sub interpreter which isn't the first interpreter created when Python is initialised, or modules which are not compatible with Python sub interpreters being destroyed and the module then being used in a new Python sub interpreter. https://code.google.com/p/modwsgi/wiki/ConfigurationDirectives#WSGIApplicationGroup"

Use the WSGIApplicationGroup directive within the virtualhost definition i.e WSGIApplicationGroup %{GLOBAL}



