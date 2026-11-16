from apcs.report import plate_html
from apcs.primersets import load_config


def test_low_depth_picks_appear_in_headline_and_untrusted_text_is_escaped():
    well = dict(well='a1', verdict='PICK (low depth)', n_reads=10, assigned_reads=12,
                n_informative_reads=10, purity=1, best_ref='</script><img src=x>',
                label='x', note='review', breadth_pct=95)
    html = plate_html(load_config(), {'a1': well}, {'title': '<script>bad()</script>',
                                                 'footer': '<img src=x>'})
    assert 'Pick <b>a1</b>' in html
    assert '<script>bad()' not in html and '</script><img' not in html
    assert 'data-well="a1"' in html and 'Aureka' in html
