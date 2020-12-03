FROM sagemath/sagemath:9.1

# Inspired from https://mybinder.readthedocs.io/en/latest/dockerfile.html#preparing-your-dockerfile
# Go to: https://nbviewer.jupyter.org/github/Thierry-Dumont/RKkit/
# then,click on "Execute on Binder" (upper right corner).

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
