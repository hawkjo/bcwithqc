"""Fixtures for the config builder's browser tests.

The tests drive gui/barcode_configurator.html in Chromium, Firefox and WebKit
through pytest-playwright (the "gui-test" dependency group). Without
Playwright this folder is skipped, so bcwithqc's own test run is unaffected.

Every page is a fresh copy of the file opened from disk with all network
requests blocked. A test fails if the page logs a JavaScript error, unless it
is marked allow_console_errors. Tests wait for data-busy="false" on #app
instead of sleeping.
"""

import importlib.util
import shutil

import pytest
from gui_support import GUI_PATH

PLAYWRIGHT_AVAILABLE = all(
    importlib.util.find_spec(name) is not None
    for name in ("playwright", "pytest_playwright")
)
if not PLAYWRIGHT_AVAILABLE:
    collect_ignore_glob = ["test_*.py"]

ENGINES = ["chromium", "firefox", "webkit"]
READY = "document.getElementById('app') && document.getElementById('app').dataset.busy === 'false'"


def pytest_configure(config):
    """Run every engine unless --browser picks some, and register this folder's markers."""
    if PLAYWRIGHT_AVAILABLE and not getattr(config.option, "browser", ["unset"]):
        config.option.browser = list(ENGINES)
    config.addinivalue_line(
        "markers",
        "allow_console_errors: the test expects the page to log JavaScript errors",
    )
    config.addinivalue_line("markers", "bcwithqc: tier F, which runs bcwithqc itself")


def pytest_report_header(config):
    """Say why the GUI tests did not run, when they can't."""
    if not PLAYWRIGHT_AVAILABLE:
        return "tests/gui: Playwright is not installed, so the GUI tests are skipped (pip install --group gui-test)."
    return None


@pytest.fixture(scope="session")
def gui_url(tmp_path_factory):
    """A file:// URL of a private copy of the page, so no test can touch the original."""
    copy = tmp_path_factory.mktemp("gui") / GUI_PATH.name
    shutil.copyfile(GUI_PATH, copy)
    return copy.as_uri()


def block_network(route):
    """Let file:// loads through and abort everything else."""
    if route.request.url.startswith("file:"):
        route.continue_()
    else:
        route.abort()


class PageErrors:
    """JavaScript errors and failed requests a page has logged."""

    def __init__(self, page):
        self.messages = []
        page.on(
            "console",
            lambda message: (
                self.messages.append(message.text) if message.type == "error" else None
            ),
        )
        page.on("pageerror", lambda error: self.messages.append(str(error)))
        page.on(
            "requestfailed",
            lambda request: self.messages.append("request blocked: " + request.url),
        )


def open_app(browser, url, **context_args):
    """Open the page in a new context with the network blocked; return (context, page, errors)."""
    context = browser.new_context(accept_downloads=True, **context_args)
    context.route("**/*", block_network)
    page = context.new_page()
    errors = PageErrors(page)
    page.goto(url)
    page.wait_for_function(READY)
    return context, page, errors


@pytest.fixture(scope="session")
def shared_app(browser, gui_url):
    """One page per engine for the pure-function tests, which don't change page state."""
    context, page, errors = open_app(browser, gui_url)
    assert not errors.messages, (
        f"JavaScript errors while loading the page: {errors.messages}"
    )
    yield page, errors
    context.close()


class Core:
    """Calls into window.barcode_core on a page."""

    def __init__(self, page):
        self.page = page

    def call(self, name, *args):
        """Call barcode_core[name](*args) and return its JSON-serialisable result (awaiting promises)."""
        return self.page.evaluate(
            "([name, args]) => window.barcode_core[name](...args)", [name, list(args)]
        )

    def run(self, script, arg=None):
        """Evaluate a JavaScript function "(core, arg) => ..." with window.barcode_core and arg."""
        return self.page.evaluate(
            "(arg) => (" + script + ")(window.barcode_core, arg)", arg
        )


@pytest.fixture
def core(shared_app, request):
    """The core functions in the shared page; fails the test if it causes a JavaScript error."""
    page, errors = shared_app
    before = len(errors.messages)
    yield Core(page)
    new_errors = errors.messages[before:]
    if request.node.get_closest_marker("allow_console_errors") is None:
        assert not new_errors, f"JavaScript errors: {new_errors}"


@pytest.fixture
def app(browser, gui_url, request):
    """A fresh page for UI tests; fails the test if the page logs a JavaScript error."""
    context, page, errors = open_app(browser, gui_url)
    yield page
    context.close()
    if request.node.get_closest_marker("allow_console_errors") is None:
        assert not errors.messages, f"JavaScript errors: {errors.messages}"
