def disableStackLimit():
    try:
        import resource
        resource.setrlimit(resource.RLIMIT_STACK, [resource.RLIM_INFINITY, resource.RLIM_INFINITY])
    except ImportError:
        pass # not a unix machine
    except ValueError:
        # could not increase stack limit
        __import__("warnings").warn("Could not change stack limit. Let's hope for the best.")
