"""
common -- plumbing every window shares.

``jobs`` runs one long operation off the GUI thread, ``notify`` shows
messages, ``gcguard`` keeps garbage collection out of Qt calls, and
``screen`` lets a tall window fit the screen it opens on, and
``platform_style`` evens out what differs between macOS and Windows.
"""
