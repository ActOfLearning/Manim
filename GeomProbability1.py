from manimlib.imports import *
import numpy as np

def Circle_Inversion(obj, about_point = ORIGIN, inversion_radius = 1.0):
    def inverse_of_a_point(pt0):
        norm_of_pt = np.dot(pt0, pt0)
        if norm_of_pt == 0:
            norm_of_pt += 0.0001
        return pt0 / norm_of_pt
    def inverse_of_a_circle(obj):
        center_of_the_circle = obj.get_center()
        radius_of_the_circle = np.linalg.norm(center_of_the_circle - obj.get_start())
        inversion_factor = np.dot(center_of_the_circle, center_of_the_circle) - radius_of_the_circle ** 2
        if inversion_factor == 0:
            res_object = Line(10 * UP, 10 * DOWN)
            res_object.move_to(RIGHT * (1 / 2 / radius_of_the_circle))
            rot_angle = np.angle(center_of_the_circle[0] + center_of_the_circle[1] * 1j)
            res_object.rotate_about_origin(rot_angle)
            return res_object
        inversion_factor = 1.0 / inversion_factor
        return Circle(radius = abs(inversion_factor) * radius_of_the_circle).move_to(inversion_factor * center_of_the_circle).set_color(TEAL)
    def circum_center(pt0):
        pt = pt0[:]
        side_a, side_b, side_c = pt[2] - pt[1], pt[0] - pt[2], pt[1] - pt[0]
        side_a, side_b, side_c = np.dot(side_a, side_a), np.dot(side_b, side_b), np.dot(side_c, side_c)
        circum_x, circum_y, circum_z = side_a * (side_b + side_c - side_a), side_b * (side_c + side_a - side_b), side_c * (side_a + side_b - side_c)
        return (circum_x * pt[0] + circum_y * pt[1] + circum_z * pt[2]) / (circum_x + circum_y + circum_z)
    def inverted_arc(c0, r0, ang0):
        return Arc(ang0[0], ang0[1] - ang0[0] + (2 * PI if ang0[1] < ang0[0] else 0), radius = r0, arc_center = c0)
    def arg_of_a_point(pt0):
        return np.angle(pt0[0] + 1j * pt0[1])
    def inverse_of_a_line(obj):
        line_st, line_en = obj.get_start(), obj.get_end()
        inv_line_start = inverse_of_a_point(line_st)
        inv_line_end = inverse_of_a_point(line_en)
        if abs(np.linalg.det([line_st, line_en, OUT])) >= 0.001:
            circum = circum_center([ORIGIN, inv_line_start, inv_line_end])
            circum_norm = np.linalg.norm(circum)
            temp_mat = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
            temp_mat = np.dot(np.dot(temp_mat, np.transpose(line_st)), line_en)
            return inverted_arc(circum, circum_norm, [arg_of_a_point(inv_line_start - circum), arg_of_a_point(inv_line_end - circum)] if temp_mat > 0 else [arg_of_a_point(inv_line_end - circum), arg_of_a_point(inv_line_start - circum)])
        else:
            if np.dot(line_st, line_en) > 0:
                return Line(inv_line_start, inv_line_end)
            else:
                return VGroup(Line(inv_line_start, 6 * inv_line_start / (0.1 + np.linalg.norm(inv_line_start))), Line(inv_line_end, 6 * inv_line_end / (0.1 + np.linalg.norm(inv_line_end))))
    def inverse_of_a_polygon(obj):
        polygon_vertices = obj.get_points()
        inverted_points = VGroup()
        for i in range(1, len(polygon_vertices)):
            inverted_points.add(inverse_of_a_line(Line(polygon_vertices[i], polygon_vertices[i - 1])))
        return inverted_points.set_color(TEAL)
    shifted_obj = obj.copy()
    shifted_obj.shift(-about_point)
    shifted_obj.scale(1 / inversion_radius, about_point = ORIGIN)
    if type(obj).__name__ == "Circle":
        res_obj = inverse_of_a_circle(shifted_obj)
    elif type(obj).__name__ == "Dot":
        res_obj = Dot(radius = DEFAULT_DOT_RADIUS / inversion_radius).move_to(inverse_of_a_point(shifted_obj.get_center()))
    elif type(obj).__name__ == "VGroup":
        res_obj = VGroup()
        for objl in shifted_obj:
            res_obj.add(Circle_Inversion(objl))
    else:
        res_obj = inverse_of_a_polygon(shifted_obj)
    res_obj.scale(inversion_radius, about_point = ORIGIN)
    res_obj.shift(about_point)
    return res_obj.set_color(TEAL)

def ParallalLine(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            return Line(p_cent, p_cent + normalized_line)
    return res

def PerpendicularFoot(pt0, line0):
    res = VGroup()
    if type(pt0).__name__ == "Dot" and type(line0).__name__ == "Line":
        line_start, line_end = line0.get_start_and_end()
        p_cent = pt0.get_center()
        if np.linalg.norm(line_end - line_start) > 0:
            normalized_line = (line_end - line_start) / np.linalg.norm(line_end - line_start)
            line_start += np.dot(p_cent - line_start, normalized_line) * normalized_line
            return Dot().move_to(line_start)
    return res

def PerpendicularLine(pt0, line0):
    line_start, line_end = line0.get_start_and_end()
    line_unit = (line_end - line_start) / np.linalg.norm(line_end - line_start)
    perp_unit = np.cross(line_unit, OUT)
    perp_foot = PerpendicularFoot(pt0, line0).get_center()
    return Line(perp_foot, perp_unit + perp_foot)

def DottoCircle(obj1):
    cen1, rad1 = ORIGIN, 0.0
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1 = obj1.get_center()
        if type(obj1).__name__ == "Dot":
            rad1 = 0.0
        else:
            rad1 = np.linalg.norm(cen1 - obj1.get_start())
    return cen1, rad1

def RadicalAxis(obj1, obj2):
    if type(obj1).__name__ == "Dot" or type(obj1).__name__ == "Circle":
        cen1, rad1 = DottoCircle(obj1)
    if type(obj2).__name__ == "Dot" or type(obj2).__name__ == "Circle":
        cen2, rad2 = DottoCircle(obj2)
    dist_cent = np.linalg.norm(cen1 - cen2)
    if dist_cent > 0:
        res = rad1 ** 2 - rad2 ** 2
        res = dist_cent + res / dist_cent
        res /= 2
        res /= dist_cent
        res = cen1 + res * (cen2 - cen1)
        return PerpendicularLine(Dot().move_to(res), Line(cen1, cen2))

def LineSegIntersection(obj1, obj2):
    res_obj = VGroup()
    obj1_start, obj1_end = obj1.get_start_and_end()
    a = obj1_end - obj1_start
    obj2_start, obj2_end = obj2.get_start_and_end()
    b = obj2_end - obj2_start
    c = obj2_start - obj1_start
    acrossb = np.cross(a, b)
    bcrossa = np.cross(b, a)
    acrossc = np.cross(a, c)
    ccrossb = np.cross(c, b)
    acrossb_dot = np.dot(acrossb, acrossb)
    bcrossa_dot = np.dot(bcrossa, bcrossa)
    if acrossb_dot > 0:
        s = np.dot(ccrossb, acrossb) / acrossb_dot
        t = np.dot(acrossc, bcrossa) / bcrossa_dot
        if 0 <= s <= 1 and 0 <= t <= 1:
            res_obj.add(Dot(obj1_start + s * a))
    return res_obj

def RayIntersection(obj1, obj2):
    #needs more work
    res_obj = VGroup()
    obj1_start, obj1_end = obj1.get_start_and_end()
    ab = obj1_end - obj1_start
    ab_unit = ab / np.linalg.norm(ab)
    obj2_start, obj2_end = obj2.get_start_and_end()
    cd = obj2_end - obj2_start
    cd_unit = cd / np.linalg.norm(cd)
    ac = obj2_start - obj1_start
    v = np.cross(ab_unit, cd_unit)
    if np.linalg.norm(v) > 0:
        v_unit = v / np.linalg.norm(v)
        v_sq = np.dot(v, v)
        s1 = np.dot(ac, np.cross(cd_unit, v_unit)) / v_sq
        s2 = np.dot(ac, np.cross(ab_unit, v_unit)) / v_sq
        print(s1, s2, obj2_start, obj2_end)
        if 0 < s1 < 1 and 0 < s2 < 1:
            return (obj1_start + ab_unit * s1 + obj2_start + cd_unit * s2) / 2
    return res_obj

def Intersection(*args):
    res_obj = VGroup()
    obj1, obj2 = args[0], args[1]
    if type(obj1).__name__ == "Line":
        obj1_start, obj1_end = obj1.get_start_and_end()
        ab = obj1_end - obj1_start
        ab_unit = ab / np.linalg.norm(ab)
        if type(obj2).__name__ == "Line":
            obj2_start, obj2_end = obj2.get_start_and_end()
            cd = obj2_end - obj2_start
            ac = obj2_start - obj1_start
            abcd = np.cross(ab, cd)
            abcd_n = np.dot(abcd, abcd)
            if abcd_n != 0:
                abcd = abcd / abcd_n
                abcd = np.dot(np.cross(ac, cd), abcd)
                abcd = obj1_start + ab * abcd
                res_obj.add(Dot(abcd))
        elif type(obj2).__name__ == "Circle":
            cen = Dot().move_to(obj2.get_center())
            p_foot = PerpendicularFoot(cen, obj1)
            p_dist = np.linalg.norm(cen.get_center() - p_foot.get_center())
            rad = np.linalg.norm(obj2.get_center() - obj2.get_start())
            if p_dist == rad:
                res_obj.add(p_foot)
            elif p_dist < rad:
                chord_dist = rad * rad - p_dist * p_dist
                chord_dist = np.sqrt(chord_dist)
                res_obj.add(Dot().move_to(p_foot.get_center() + chord_dist * ab_unit), Dot().move_to(p_foot.get_center() - chord_dist * ab_unit))
    elif type(obj1).__name__ == "Circle":
        cen1, cen2 = obj1.get_center(), obj2.get_center()
        dist_cen = np.linalg.norm(cen1 - cen2)
        rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
        if abs(dist_cen - rad1 - rad2) < 1e-3:
            res_obj.add(Dot().move_to((rad2 * cen1 + rad1 * cen2) / (rad1 + rad2)))
        elif dist_cen < rad1 + rad2:
            p_foot = RadicalAxis(obj1, obj2)
            #p_foot = PerpendicularLine(p_foot, Line(cen1, cen2))
            return Intersection(p_foot, obj1)
    return res_obj

def Midpoint(*args):
    if type(args[0]).__name__ == "Line":
        line_start, line_end = args[0].get_start_and_end()
        mid_pt = (line_start + line_end) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        mid_pt = (args[0].get_center() + args[1].get_center()) / 2
        mid_pt = Dot().move_to(mid_pt)
        return mid_pt

def PerpendicularBisector(*args):
    mid_pt = Midpoint(*args)
    if type(args[0]).__name__ == "Line":
        print("It's a line..")
        return PerpendicularLine(mid_pt, args[0])
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot":
        connec_line = Line(args[0].get_center(), args[1].get_center())
        return PerpendicularLine(mid_pt, connec_line)

def AngularBisector(*args):
    if type(args[0]).__name__ == "Line" and type(args[1]).__name__ == "Line":
        intersec = Intersection(*args)
        if len(list(intersec)) > 0:
            return AngularBisector(Dot().move_to(args[0].get_end()), intersec, Dot().move_to(args[1].get_end()))
        else:
            return VGroup()
    elif type(args[0]).__name__ == "Dot" and type(args[1]).__name__ == "Dot" and type(args[2]).__name__ == "Dot":
        pt_a = args[0].get_center()
        pt_b = args[1].get_center()
        pt_c = args[2].get_center()
        vec_ba = pt_a - pt_b
        vec_bc = pt_c - pt_b
        return Line(pt_b, vec_ba + vec_bc)

def RadicalCenter(*args):
    if type(args[0]).__name__ == "Circle" and type(args[1]).__name__ == "Circle" and type(args[2]).__name__ == "Circle":
        radic12 = RadicalAxis(args[0], args[1])
        radic23 = RadicalAxis(args[1], args[2])
        return Intersection(radic12, radic23)
    return Dot()

def HomotheticCenter(obj1, obj2):
    cen1, cen2 = obj1.get_center(), obj2.get_center()
    rad1, rad2 = np.linalg.norm(cen1 - obj1.get_start()), np.linalg.norm(cen2 - obj2.get_start())
    res_obj = VGroup()
    res_obj.add(Dot().move_to((rad1 * cen2 + rad2 * cen1) / (rad1 + rad2)))
    if rad1 != rad2:
        res_obj.add(Dot().move_to((rad1 * cen2 - rad2 * cen1) / (rad1 - rad2)))
    return res_obj

def CircleTangentLines(obj, pt):
    cen = obj.get_center()
    rad = np.linalg.norm(cen - obj.get_start())
    dist = np.linalg.norm(cen - pt.get_center())
    res_obj = VGroup()
    if dist < rad:
        return res_obj.add(pt)
    temp_c = Circle(arc_center = (cen + pt.get_center()) / 2, radius = dist / 2)
    temp_c = Intersection(temp_c, obj)
    for i in temp_c:
        res_obj.add(Line(pt.get_center(), i.get_center()))
    return res_obj

class Angle(VGroup):

    CONFIG = {
        'radius': 0.5,
        'color': RED,
        'opacity': 0.4,
        'stroke_width': 5,
        'below_180': True,
    }

    def __init__(self, A, O, B, **kwargs):

        VMobject.__init__(self, **kwargs)
        OA, OB = A-O, B-O
        if self.below_180:
            theta = np.angle(complex(*OA[:2])/complex(*OB[:2])) # angle of OB to OA
        else:
            theta = TAU + np.angle(complex(*OA[:2])/complex(*OB[:2]))

        if abs(abs(theta) - PI / 2) <= 1e-3:
            self.add(Elbow(angle = np.angle(complex(*OB[:2])), color = self.color, width = self.radius).shift(O))
        else:
            self.add(Arc(start_angle=Line(O, B).get_angle(), angle=theta, radius=self.radius/2,
                     stroke_width=100 * self.radius, color=self.color, arc_center=O).set_stroke(opacity=self.opacity))
            self.add(Arc(start_angle=Line(O, B).get_angle(), angle=theta, radius=self.radius,
                     stroke_width=self.stroke_width, color=self.color, arc_center=O))

class ProblemIntroduction_Obsolete(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10,
        "x_max": 10,
        "x_axis_width": 20,
        "x_tick_frequency": 1,
        "x_leftmost_tick": None,  # Change if different from x_min
        "x_labeled_nums": None,
        "x_axis_label": "",
        "y_min": -5,
        "y_max": 5,
        "y_axis_height": 10,
        "y_tick_frequency": 1,
        "y_bottom_tick": None,  # Change if different from y_min
        "y_labeled_nums": None,
        "y_axis_label": "",
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
        "exclude_zero_label": True,
        "default_graph_colors": [YELLOW, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": YELLOW,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.25,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": YELLOW,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
        "given_circle_color": RED,
        "incircle_color": PURPLE,
        "excircle_color": YELLOW,
        "line_of_centers": GREY,
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        geoprob = TextMobject("\\text{\\underline{Geometric Probability}}").scale(1.125)
        self.play(Write(geoprob))
        self.play(geoprob.move_to, 0.8125 * TOP)
        tri = VMobject(stroke_color = YELLOW)
        a, b, c = 3 * UP + 2 * RIGHT, ORIGIN, 5 * RIGHT
        adot, bdot, cdot = Dot(a), Dot(b), Dot(c)
        alabel, blabel, clabel = TextMobject("A"), TextMobject("B"), TextMobject("C")
        alabel.next_to(adot, direction = UR, buff = 0.04)
        blabel.next_to(bdot, direction = DOWN, buff = 0.05)
        clabel.next_to(cdot, direction = DOWN, buff = 0.05)
        labelgrp = VGroup(alabel, blabel, clabel)
        ab = Line(a, b)
        bc = Line(b, c)
        ca = Line(c, a)
        ddot, edot, fdot = PerpendicularFoot(adot, bc), PerpendicularFoot(bdot, ca), PerpendicularFoot(cdot, ab)
        dlabel, elabel, flabel = TextMobject("D"), TextMobject("E"), TextMobject("F")
        dlabel.next_to(ddot, direction = DOWN, buff = 0.05)
        elabel.next_to(edot, direction = RIGHT, buff = 0.05)
        flabel.next_to(fdot, direction = LEFT, buff = 0.05)
        def_label = VGroup(dlabel, elabel, flabel)
        dangle, eangle, fangle = Angle(a, ddot.get_center(), c), Angle(a, edot.get_center(), b), Angle(c, fdot.get_center(), b)
        dangle.scale_in_place(0.5, about_point = ddot.get_center())
        eangle.scale_in_place(0.5, about_point = edot.get_center())
        fangle.scale_in_place(0.5, about_point = fdot.get_center())
        angle_grp = VGroup(dangle, eangle, fangle)
        ad, be, cf = Line(adot.get_center(), ddot.get_center()), Line(bdot.get_center(), edot.get_center()), Line(cdot.get_center(), fdot.get_center())
        altitude_grp = VGroup(ad, be, cf)
        altitude_grp.set_color(BLUE)
        tri.set_points_as_corners([a, b, c, a])
        self.play(
            self.camera_frame.shift, 2 * RIGHT + 2 * UP,
            geoprob.shift, 2 * RIGHT + 2 * UP,
            FadeIn(tri),
            FadeIn(labelgrp))
        self.wait()
        self.play(*[Write(obj) for obj in [ddot, edot, fdot]])
        self.wait()
        self.play(ShowCreation(altitude_grp), Write(def_label))
        self.wait()
        self.play(Write(angle_grp))
        self.wait()
        self.play(*[FadeOut(obj) for obj in [angle_grp, altitude_grp, def_label, ddot, edot, fdot]])
        self.wait()
        rdots, rlines = VGroup(), VGroup()
        def triintersector(r_line, obj_ab, obj_bc, obj_ca):
                for k in [obj_ab, obj_bc, obj_ca]:
                    tempobj = LineSegIntersection(r_line, k)
                    if len(tempobj) > 0:
                        return tempobj[0]
        def rptangle():
            ac_rand = random.random()
            ab_rand = random.random()
            if ac_rand + ab_rand > 1:
                ac_rand, ab_rand = 1 - ac_rand, 1 - ab_rand
            rpt = a + (c - a) * ac_rand + (b - a) * ab_rand
            #rdot = Dot(rpt).scale_in_place(0.5)
            rangle = random.random() * 2 * PI
            return rpt, rangle
        def dot_creator(pnt, angl):
            t, a = rptangle()
            rlin = Line(t, t + 5 * np.array([np.cos(a), np.sin(a), 0]))
            pt = triintersector(rlin, ab, bc, ca)
            #rline = Line(t, pt.get_center(), stroke_width = 0.5)
            rline = Line(t, pt.get_center())
            return rline
        m, n = rptangle()
        rline = dot_creator(m, n)
        rlines.add(rline)
        self.play(Write(rline))
        self.wait()
        for _ in range(50):
            #self.remove(rdot, rline)
            m, n = rptangle()
            rline = dot_creator(m, n)
            rlines.add(rline)
            self.play(Write(rline), run_time = 0.2)
        self.wait()
        '''for _ in range(700):
            self.remove(rline)
            m, n = rptangle()
            rline = dot_creator(m, n)
            rlines.add(rline)'''
        self.play(*[FadeOut(obj) for obj in [rlines, tri, labelgrp]])
        self.wait()
        requisites = TextMobject("\\text{\\underline{Pre - Requisites}}").scale(1.125)
        requisites.move_to(geoprob.get_center())
        self.play(Transform(geoprob, requisites))
        self.wait()
        text_scale = 0.75
        reqpts = VGroup(
            TextMobject("""$\\bullet$ Given a random variable $X$ with pdf $f(x)$, $\\mathbb{E}(g(X))=\\displaystyle\\int\\limits_{-\\infty}^\\infty g(x)f(x)\\,dx$""").scale(text_scale),
            TextMobject("""$\\circ$ If $X \\sim U(0, 2\\pi)$, then $\\mathbb{E}(g(X))=\\displaystyle\\frac{1}{2\\pi}\\int\\limits_{0}^{2\\pi} g(x)\\,dx$""").scale(text_scale),
            TextMobject("""$\\circ$ If $P(x,y) \\sim \\triangle ABC$, then $\\mathbb{E}(g(P))=\\displaystyle\\frac{1}{\\text{Area}(\\triangle ABC)}\\iint\\limits_{(x,y)\\in\\triangle ABC} g(P)\\,dydx$""").scale(text_scale),
            TextMobject("""$\\bullet$ $\\displaystyle\\int\\limits_{0}^{\\tan^{-1}(x)} \\text{sec}(\\theta) \\,d\\theta=\\sinh^{-1}(x)$""").scale(text_scale),
            TextMobject("""$\\bullet$ $\\mathbb{E}(Y)=\\mathbb{E}(\\mathbb{E}(Y|X))$ (or) $\\mathbb{E}(Y)=\\displaystyle\\int\\limits_{-\\infty}^\\infty \\mathbb{E}(Y|X=x)f(x)\\,dx$""").scale(text_scale),
        )
        reqpts.arrange(direction = DOWN, buff = 0.01)
        for obj in reqpts[1:]:
            obj.align_to(reqpts[0], LEFT)
        reqpts[1].shift(0.75 * RIGHT + 0.25 * UP)
        reqpts[2].shift(0.75 * RIGHT + 0.25 * UP)
        reqpts.shift(1.5 * RIGHT + 1.5 * UP)
        for text in reqpts:
            self.play(Write(text))
            self.wait()
        self.wait(5)

class Problemsolution(GraphScene, MovingCameraScene):
    CONFIG = {
        "x_min": -10,
        "x_max": 10,
        "x_axis_width": 20,
        "x_tick_frequency": 1,
        "x_leftmost_tick": None,  # Change if different from x_min
        "x_labeled_nums": None,
        "x_axis_label": "",
        "y_min": -5,
        "y_max": 5,
        "y_axis_height": 10,
        "y_tick_frequency": 1,
        "y_bottom_tick": None,  # Change if different from y_min
        "y_labeled_nums": None,
        "y_axis_label": "",
        "axes_color": BLACK,
        "graph_origin": ORIGIN,
        "exclude_zero_label": True,
        "default_graph_colors": [YELLOW, GREEN, YELLOW],
        "default_derivative_color": GREEN,
        "default_input_color": YELLOW,
        "default_riemann_start_color": YELLOW,
        "default_riemann_end_color": GREEN,
        "area_opacity": 0.8,
        "num_rects": 50,
        "dot_kwargs": {
            "radius": 0.25,
        },
        "line_kwargs": {
            "stroke_width": 2,
        },
        "fill_triangle_kwargs": {
            # "fill_color": YELLOW,
            "fill_opacity": .5,
            "stroke_width": 0,
        },
        "given_circle_color": RED,
        "incircle_color": PURPLE,
        "excircle_color": YELLOW,
        "line_of_centers": GREY,
    }

    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = NumberPlane(y_min = -3, y_max = 10)
        self.play(
            self.camera_frame.scale, 1,
        )
        tri = VMobject(stroke_color = YELLOW)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        adot, bdot, cdot = Dot(a), Dot(b), Dot(c)
        text_scale = 0.6875
        alabel, blabel, clabel = TextMobject("A").scale(text_scale), TextMobject("B").scale(text_scale), TextMobject("C").scale(text_scale)
        alabel.next_to(adot, direction = UR, buff = 0.04)
        blabel.next_to(bdot, direction = DR, buff = 0.05)
        clabel.next_to(cdot, direction = DL, buff = 0.05)
        labelgrp = VGroup(alabel, blabel, clabel)
        ab = Line(a, b)
        bc = Line(b, c)
        ca = Line(c, a)
        ddot, edot, fdot = PerpendicularFoot(adot, bc), PerpendicularFoot(bdot, ca), PerpendicularFoot(cdot, ab)
        dlabel, elabel, flabel = TextMobject("D").scale(text_scale), TextMobject("E").scale(text_scale), TextMobject("F").scale(text_scale)
        dlabel.next_to(ddot, direction = DOWN, buff = 0.05)
        elabel.next_to(edot, direction = RIGHT, buff = 0.05)
        flabel.next_to(fdot, direction = LEFT, buff = 0.05)
        def_label = VGroup(dlabel, elabel, flabel)
        dangle, eangle, fangle = Angle(a, ddot.get_center(), c), Angle(a, edot.get_center(), b), Angle(c, fdot.get_center(), b)
        dangle.scale_in_place(0.5, about_point = ddot.get_center())
        eangle.scale_in_place(0.5, about_point = edot.get_center())
        fangle.scale_in_place(0.5, about_point = fdot.get_center())
        angle_grp = VGroup(dangle, eangle, fangle)
        ad, be, cf = Line(adot.get_center(), ddot.get_center()), Line(bdot.get_center(), edot.get_center()), Line(cdot.get_center(), fdot.get_center())
        altitude_grp = VGroup(ad, be, cf)
        altitude_grp.set_color(BLUE)
        tri.set_points_as_corners([a, b, c, a])
        def triintersector(r_line, obj_ab, obj_bc, obj_ca):
                for k in [obj_ab, obj_bc, obj_ca]:
                    tempobj = LineSegIntersection(r_line, k)
                    if len(tempobj) > 0:
                        return tempobj[0]
        def rptangle(obj_a, obj_b, obj_c):
            ac_rand = random.random()
            ab_rand = random.random()
            if ac_rand + ab_rand > 1:
                ac_rand, ab_rand = 1 - ac_rand, 1 - ab_rand
            rpt = obj_a + (obj_c - obj_a) * ac_rand + (obj_b - obj_a) * ab_rand
            #rdot = Dot(rpt).scale_in_place(0.5)
            rangle = random.random() * 2 * PI
            return rpt, rangle
        def dot_creator(pnt, angl, obj_ab, obj_bc, obj_ca):
            rlin = Line(pnt, pnt + 5 * np.array([np.cos(angl), np.sin(angl), 0]))
            pt = triintersector(rlin, obj_ab, obj_bc, obj_ca)
            if pt == None:
                return Line(ORIGIN, ORIGIN + RIGHT / 1000)
            rline = Line(pnt, pt.get_center())
            return rline
        self.play(
            self.camera_frame.shift, 2.5 * UP,
            FadeIn(tri),
            FadeIn(labelgrp),
        )
        self.wait()
        self.play(Write(nplane))
        self.wait()
        self.play(
            Write(altitude_grp),
            Write(def_label)
        )
        hlengths = VGroup(
            TextMobject("$h_a=|AD|$").scale(text_scale),
            TextMobject("$h_b=|BE|$").scale(text_scale),
            TextMobject("$h_c=|CF|$").scale(text_scale),
        )
        hlengths.arrange(direction = DOWN, buff = 0.1)
        hlengths.move_to(5 * RIGHT + 3.5 * UP)
        hrect = BackgroundRectangle(hlengths)
        hset = VGroup(hrect, hlengths)
        self.add_foreground_mobject(hset)
        self.play(
            Write(hset),
        )
        self.wait()
        self.play(
            FadeOut(altitude_grp),
            FadeOut(def_label),
            FadeOut(hset)
        )
        self.wait()
        m, n = 1.75 * RIGHT + 1.25 * UP, (PI / 4) * 1.01342687
        rline = dot_creator(m, n, ab, bc, ca)
        llabel = TextMobject("L").scale(text_scale)
        llabel.next_to(rline, UP, buff = -0.5)
        self.play(Write(rline), Write(llabel))
        self.wait()
        eqn_scale = 0.75
        leqn = VGroup(
            TextMobject("$L=g(x,y,t)$").scale(eqn_scale),
            TextMobject("$\\mathbb{E}(L)=\\displaystyle \\int\\limits_{0}^{2\\pi}\\int\\int g(x,y,t)\\,\\frac{dydx}{\\triangle ABC}\\frac{dt}{2\\pi}$").scale(eqn_scale),
        )
        leqn.arrange(direction = DOWN, buff = 0.5)
        leqn.move_to(5 * LEFT + 3 * UP)
        for obj in leqn[1:]:
            obj.align_to(leqn[0], LEFT)
        ddot = Dot().scale(text_scale).move_to(m)
        xylabel = TextMobject("$(x,y)$").scale(text_scale)
        tlabel = TexMobject("t").scale(text_scale)
        xylabel.move_to(m)
        xylabel.shift(0.5 * LEFT)
        dline = DashedLine(m, m + RIGHT).set_color(GREY)
        tang = Angle(dline.get_end(), m, rline.get_end())
        tlabel.move_to(2.5 * RIGHT + 1.5 * UP)
        self.play(
            *[Write(obj) for obj in [ddot, xylabel, tlabel, dline, tang, leqn[0]]]
        )
        self.wait()
        self.play(Write(leqn[1]))
        self.wait()
        cansolvetxt = TextMobject("Can this be solved \\\\ in closed form?")
        tbubble = ThoughtBubble()
        tbubble.pin_to(leqn[1])
        tbubble.add_content(cansolvetxt)
        tbubble.resize_to_content()
        bub_grp = VGroup(tbubble, cansolvetxt)
        self.play(Write(bub_grp))
        self.wait()
        self.play(FadeOut(bub_grp), FadeOutAndShift(leqn, direction = DOWN))
        self.wait()
        condtxt = TextMobject("Conditioning is the soul of Statistics", tex_to_color_map = {"Conditioning": YELLOW})
        condtxt.move_to(5.5 * UP)
        self.play(Write(condtxt))
        self.wait()
        self.play(FadeOut(condtxt))
        self.wait()
        setuptxt = TextMobject("Setup S1: Choose a point randomly in $\\triangle ABC$ and draw a line perpendiuclar to base $BC$ from that point.", tex_to_color_map = {"Setup S1": YELLOW}).scale(eqn_scale)
        setupbox = BackgroundRectangle(setuptxt)
        setupgrp = VGroup(setupbox, setuptxt)
        setupgrp.move_to(5.5 * UP)
        #self.add_foreground_mobject(setupgrp)
        self.play(
            Write(setupgrp),
            *[FadeOut(obj) for obj in [ddot, xylabel, tlabel, dline, tang, llabel, rline]]
        )
        self.wait()
        rlines = VGroup()
        rline = dot_creator(2.25 * UR, -PI / 2, ab, bc, ca)
        rlines.add(rline)
        self.play(Write(rline))
        #self.remove(xylabel)
        xylabel.next_to(rline.get_start(), direction = UP, buff = 0.05)
        self.play(Write(xylabel))
        self.wait()
        expeceqn = TexMobject("\\mathbb{E}_{ABC}(L|S_1)").scale(eqn_scale)
        eqnlist = VGroup(
            TexMobject("=\\displaystyle\\frac{1}{\\text{Area}(\\triangle)}\\int\\limits_{(x,y)\\in \\triangle}y \\,dydx").scale(eqn_scale),
            TexMobject("=\\displaystyle\\frac{|AD|}{3}").scale(eqn_scale),
        )
        eqnlist.arrange(DOWN)
        eqnlist.move_to(2.25 * LEFT + 3 * UP)
        expeceqn.next_to(eqnlist[0], direction = LEFT)
        expeceqn.shift(0.125 * UP)
        for obj in eqnlist[1:]:
            obj.align_to(eqnlist[0], LEFT)
        self.play(Write(expeceqn), Write(eqnlist[0]))
        self.wait()
        self.play(FadeOut(xylabel))
        '''for _ in range(50):
            m, n = rptangle(a, b, c)
            rline = dot_creator(m, -PI / 2, ab, bc, ca)
            rlines.add(rline)
            self.play(Write(rline), run_time = 0.15)
        self.wait()'''
        self.play(Write(dlabel), Write(ad), FadeOut(rlines))
        self.wait()
        self.play(Write(eqnlist[1]))
        self.wait()
        ddot = PerpendicularFoot(adot, bc).scale_in_place(text_scale)
        self.play(*[FadeOut(obj) for obj in [ad, eqnlist, expeceqn]], FadeIn(ddot))
        self.wait()
        ##########################################################################
        self.play(
            Transform(setuptxt,
            TextMobject("Setup S2: Given $M$ between $D$ and $C$, choose a point randomly in $\\triangle ABC$ and draw a line parallel to $AM$ $towards$ base $BC$.", tex_to_color_map = {"Setup S2": YELLOW, "$towards$": PINK}).scale(eqn_scale).move_to(5.5 * UP)
            )
        )
        self.wait()
        mlabel = TextMobject("M").scale(text_scale)
        mdot = Dot(3 * RIGHT).scale_in_place(text_scale)
        mlabel.next_to(mdot, direction = DOWN, buff = 0.05)
        am = DashedLine(a, mdot.get_center()).set_color(GREEN)
        self.play(*[Write(obj) for obj in [mlabel, mdot, am]])
        self.wait()
        m, n = rptangle(a, b, c)
        n = am.get_angle()
        rlines = VGroup()
        rline = dot_creator(m, n, ab, bc, ca)
        rlines.add(rline)
        self.play(Write(rline))
        plabel, qlabel = TexMobject("P(x,y)").scale(text_scale), TexMobject("Q").scale(text_scale)
        plabel.next_to(rline.get_start(), direction = UP, buff = 0.05)
        qlabel.next_to(rline.get_end(), direction = DOWN, buff = 0.05)
        self.play(Write(plabel), Write(qlabel))
        self.wait()
        xylabel.next_to(rline.get_start(), direction = UP)
        madangle = Angle(ddot.get_center(), a, mdot.get_center(), radius = 1.25)
        tlabel.move_to(1.75 * RIGHT + 2.5 * UP)
        self.play(Write(madangle), Write(tlabel), Write(ad))
        pxeqn = TextMobject("$|PQ|=\\displaystyle\\frac{y}{\\cos t}$").scale(text_scale)
        pxeqn.move_to(4.5 * LEFT + 4 * UP)
        self.play(Write(pxeqn))
        expeceqn = TexMobject("\\mathbb{E}_{ABC}(L|S_2, M)").scale(text_scale)
        eqnlist = VGroup(
            TexMobject("=\\displaystyle\\int\\limits_{(x,y)\\in \\triangle}\\frac{y}{\\cos t} \\,\\frac{dydx}{\\text{Area}(\\triangle)}").scale(text_scale),
            TexMobject("=\\displaystyle\\frac{|AD|}{3\\cos t}=\\frac{|AM|}{3}").scale(text_scale),
        )
        eqnlist.arrange(DOWN)
        eqnlist.move_to(2.25 * LEFT + 2 * UP)
        expeceqn.next_to(eqnlist[0], direction = LEFT)
        expeceqn.shift(0.125 * UP)
        for obj in eqnlist[1:]:
            obj.align_to(eqnlist[0], LEFT)
        self.play(Write(expeceqn), Write(eqnlist[0]))
        self.wait()
        self.play(Write(eqnlist[1]))
        self.wait()
        self.play(*[FadeOut(obj) for obj in [ad, expeceqn, pxeqn, eqnlist, tlabel, madangle, plabel, qlabel]])
        self.wait()
        '''for _ in range(50):
            m, n = rptangle(a, b, c)
            n = am.get_angle()
            rline = dot_creator(m, n, ab, bc, ca)
            rlines.add(rline)
            self.play(Write(rline), run_time = 2 / 50)
        self.wait()'''
        self.play(FadeOut(rlines))
        self.wait()
        ##########################################################################
        self.play(
            Transform(setuptxt,
            TextMobject("Setup S3: Given $M$ between $D$ and $C$, choose a point randomly in $\\triangle ABC$ and draw a line parallel to $AM$ $away$ from base $BC$.", tex_to_color_map = {"Setup S3": YELLOW, "$away$": PINK}).scale(eqn_scale).move_to(5.5 * UP)
            )
        )
        self.wait()
        m = 3.25 * RIGHT + 1.5 * UP
        n = am.get_angle() + PI
        rlines = VGroup()
        rline = dot_creator(m, n, ab, bc, ca)
        rlines.add(rline)
        self.play(Write(rline))
        self.wait()
        '''for _ in range(50):
            m, n = rptangle(a, b, c)
            n = am.get_angle() + PI
            rline = dot_creator(m, n, ab, bc, ca)
            rlines.add(rline)
            self.play(Write(rline), run_time = 4 / 75)
        self.wait()'''
        self.play(FadeOut(rlines))
        self.wait()
        rlines = VGroup()
        m, n = rptangle(a, mdot.get_center(), c)
        n = am.get_angle() + PI
        rline = dot_creator(m, n, ca, Line(a, mdot.get_center()), Line(mdot.get_center(), c))
        self.wait()
        pxeqn = TextMobject("$\\mathbb{E}_{ABC}(L|S_3,M)=\\displaystyle\\frac{\\triangle MAB}{\\triangle ABC}\\mathbb{E}(P \\in \\triangle MAB)+$", "$\\displaystyle \\frac{\\triangle MCA}{\\triangle ABC}\\mathbb{E}(P \\in \\triangle MCA)$").scale(0.625)
        pxeqn.arrange(DOWN, buff = 0.25)
        pxeqn[1].shift(1.75 * RIGHT)
        pxeqn.move_to(6.5 * LEFT + 4.5 * UP + pxeqn.get_center() - pxeqn.get_critical_point(UL))
        self.play(Write(pxeqn))
        #expeceqn = TexMobject("\\mathbb{E}(P \\in \\triangle AMC)=").scale(text_scale)
        eqnlist = VGroup(
            TexMobject("\\mathbb{E}(L|P \\in \\triangle MCA)=").scale(0.625),
            TexMobject("\\mathbb{E}_{MCA}(L|S_2,A)=\\displaystyle\\frac{|AM|}{3}").scale(0.625),
        )
        eqnlist[1].next_to(eqnlist[0], RIGHT)
        eqnlist[1].shift(0.0625 * UP)
        eqnlist.move_to(3.5 * LEFT + 2 * UP)
        self.wait()
        self.play(Write(eqnlist[0]))
        self.wait()
        '''rlines = VGroup()
        for _ in range(50):
            m, n = rptangle(mdot.get_center(), c, a)
            n = am.get_angle() + PI
            rline = dot_creator(m, n, Line(mdot.get_center(), c), ca, Line(a, mdot.get_center()))
            rlines.add(rline)
            self.play(Write(rline), run_time = 4 / 75)
        self.wait()'''
        self.play(Write(eqnlist[1]))
        self.wait()
        self.play(FadeOut(rlines))
        self.wait()
        '''rlines = VGroup()
        for _ in range(50):
            m, n = rptangle(mdot.get_center(), a, b)
            n = am.get_angle() + PI
            rline = dot_creator(m, n, Line(mdot.get_center(), a), ab, Line(b, mdot.get_center()))
            rlines.add(rline)
            self.play(Write(rline), run_time = 4 / 75)
        self.wait()'''
        eqnlist1 = TexMobject("\\mathbb{E}(L|P \\in \\triangle MAB)=\\mathbb{E}_{MAB}(L|S_2,A)=\\displaystyle\\frac{|AM|}{3}").scale(0.625)
        eqnlist1.align_to(eqnlist[0], LEFT)
        eqnlist1.shift(UP)
        self.play(Write(eqnlist1))
        self.wait()
        self.play(FadeOut(rlines))
        self.wait()
        restxt = TexMobject("=\\displaystyle\\frac{|AM|}{3}").scale(0.625)
        restxt.move_to(4 * LEFT + 2.5 * UP + restxt.get_center() - restxt.get_critical_point(UL))
        self.play(
            *[FadeOutAndShift(obj, direction = DOWN) for obj in [eqnlist, eqnlist1]],
            FadeIn(restxt)
        )
        self.play(*[FadeOut(obj) for obj in [restxt, pxeqn]])
        self.wait()
        ##########################################################################
        self.play(
            Transform(setuptxt,
            TextMobject("Setup S4: Given $M$ between $D$ and $C$, choose a point randomly in $\\triangle ABC$ and draw a line parallel to $AM$.", tex_to_color_map = {"Setup S4": YELLOW}).scale(eqn_scale).move_to(5.5 * UP)
            )
        )
        self.wait()
        '''rlines = VGroup()
        for _ in range(50):
            m, n = rptangle(a, b, c)
            n = am.get_angle() + (PI if random.random() <= 0.5 else 0)
            rline = dot_creator(m, n, ab, bc, ca)
            rlines.add(rline)
            self.play(Write(rline), run_time = 5 / 50)
        self.wait()'''
        assumptxt = TextMobject("Let $\\mathbb{P}(\\text{towards }BC)=p$ and ", "$\\mathbb{P}(\\text{away from }BC)=1-p$").scale(0.625)
        assumptxt.move_to(6.5 * LEFT + 4.5 * UP + assumptxt.get_center() - assumptxt.get_critical_point(UL))
        self.play(Write(assumptxt))
        lefteqn = TexMobject(
            "\\mathbb{E}(L|S_4,M)"
        ).scale(0.625)
        lefteqn.align_to(assumptxt, LEFT)
        lefteqn.shift(3.25 * UP)
        pxeqn = VGroup(
            TexMobject("=p\\mathbb{E}(\\text{towards }BC)+(1-p)\\mathbb{E}(\\text{away from }BC)").scale(0.625),
            TexMobject("=\\displaystyle p\\frac{|AM|}{3}+(1-p)\\frac{|AM|}{3}").scale(0.625),
            TexMobject("=\\displaystyle \\frac{|AM|}{3}").scale(0.625),
            TexMobject("=\\displaystyle \\frac{|AD|}{3\\cos t}=\\frac{h_a}{3\\cos t}")
        )
        pxeqn.arrange(DOWN)
        for obj in pxeqn[1:]:
            obj.align_to(pxeqn[0], LEFT)
        pxeqn.next_to(lefteqn, direction = RIGHT, aligned_edge = UL, buff = 0.9375)
        self.play(Write(lefteqn), Write(pxeqn[:3]))
        self.wait()
        self.play(
            Write(pxeqn[-1]),
            *[FadeIn(obj) for obj in [tlabel, madangle, ad]]
        )
        self.play(
            *[FadeOutAndShiftDown(obj) for obj in [rlines, pxeqn, lefteqn, assumptxt]],
            *[FadeOut(obj) for obj in [tlabel, madangle, ad, mdot, mlabel, am]]
        )
        self.wait()
        ##########################################################################
        self.play(
            Transform(setuptxt,
            TextMobject("Setup S5: Choose a point randomly in $\\triangle ABC$ and draw a line (towards or away from $BC$) such that the angle $t$ from $AD$ is $U(0,\\tan^{-1}(CD/AD))$.", tex_to_color_map = {"Setup S5": YELLOW}).scale(eqn_scale).move_to(5.5 * UP)
            )
        )
        self.wait()
        '''rlines = VGroup()
        for _ in range(50):
            m, n = rptangle(a, b, c)
            n = (-PI / 2 + random.random() * (PI + ca.get_angle() + PI / 2)) + (PI if random.random() <= 0.5 else 0)
            rline = dot_creator(m, n, ab, bc, ca)
            rlines.add(rline)
            self.play(Write(rline), run_time = 5 / 50)
        self.wait()'''
        lefteqn = TexMobject(
            "\\mathbb{E}(L|S_5,CD)"
        ).scale(0.625)
        lefteqn.move_to(4 * UP + 6.5 * LEFT + lefteqn.get_center() - lefteqn.get_critical_point(UL))
        pxeqn = VGroup(
            TexMobject("=\\mathbb{E}(\\mathbb{E}(L|S_4,t))").scale(0.625),
            TexMobject("=\\displaystyle \\int\\limits_{0}^{\\tan^{-1}(CD/AD)}\\frac{h_a}{3\\cos t}\\,\\frac{dt}{\\tan^{-1}(CD/AD)}").scale(0.625),
            TexMobject("=\\displaystyle \\frac{h_a}{3\\tan^{-1}(CD/AD)}\\text{sinh}^{-1}(CD/AD)").scale(0.625),
        )
        pxeqn.arrange(DOWN)
        for obj in pxeqn[1:]:
            obj.align_to(pxeqn[0], LEFT)
        pxeqn.next_to(lefteqn, direction = RIGHT, aligned_edge = UL, buff = 0.9375)
        self.play(Write(lefteqn))
        self.wait()
        for obj in pxeqn:
            self.play(Write(obj))
            self.wait()
        self.play(
            *[FadeOutAndShiftDown(obj) for obj in [lefteqn, pxeqn]]
        )
        self.wait()
        ##########################################################################
        self.play(
            Transform(setuptxt,
            TextMobject("Setup S6: Choose a point randomly in $\\triangle ABC$ and draw a line (towards or away from $BC$) such that it's angle $t$ from $AC$ is $U(0,\\angle A)$.", tex_to_color_map = {"Setup S6": YELLOW}).scale(eqn_scale).move_to(5.5 * UP)
            )
        )
        self.wait()
        '''rlines = VGroup()
        for _ in range(50):
            m, n = rptangle(a, b, c)
            n = (ab.get_angle() + random.random() * (PI + ca.get_angle() - ab.get_angle())) + (PI if random.random() <= 0.5 else 0)
            rline = dot_creator(m, n, ab, bc, ca)
            rlines.add(rline)
            self.play(Write(rline), run_time = 5 / 50)
        self.wait()'''
        lefteqn = TexMobject(
            "\\mathbb{E}(L|S_6,A)"
        ).scale(0.625)
        lefteqn.move_to(4 * UP + 6.5 * LEFT + lefteqn.get_center() - lefteqn.get_critical_point(UL))
        pxeqn = VGroup(
            TexMobject("=\\displaystyle \\frac{\\tan^{-1}\\left(\\frac{BD}{AD}\\right)}{\\angle A}\\mathbb{E}(L|S_5,BD)+\\frac{\\tan^{-1}\\left(\\frac{CD}{AD}\\right)}{\\angle A}\\mathbb{E}(L|S_5,CD)").scale(0.625),
            TexMobject("=\\displaystyle \\frac{h_a}{3\\angle A}\\left(\\text{sinh}^{-1}\\left(\\frac{BD}{AD}\\right)+\\text{sinh}^{-1}\\left(\\frac{CD}{AD}\\right)\\right)").scale(0.625),
            TexMobject("=\\displaystyle \\frac{h_a}{3\\angle A} \\text{sinh}^{-1}\\left(\\frac{a}{p}\\frac{p-a}{s-a}\\right)").scale(0.625),
        )
        pxeqn.arrange(DOWN)
        for obj in pxeqn[1:]:
            obj.align_to(pxeqn[0], LEFT)
        pxeqn.next_to(lefteqn, direction = RIGHT, aligned_edge = UL, buff = 0.9375)
        self.play(Write(lefteqn))
        self.wait()
        for obj in pxeqn:
            self.play(Write(obj))
            self.wait()
        self.play(
            *[FadeOutAndShiftDown(obj) for obj in [lefteqn, pxeqn]]
        )
        self.wait()
        ##########################################################################
        self.play(
            Transform(setuptxt,
            TextMobject("Setup S7: Choose a point randomly in $\\triangle ABC$ and draw a line at a random angle from that point until it meets the triangle.", tex_to_color_map = {"Setup S7": YELLOW}).scale(eqn_scale).move_to(5.5 * UP)
            )
        )
        self.wait()
        '''rlines = VGroup()
        for _ in range(50):
            m, n = rptangle(a, b, c)
            rline = dot_creator(m, n, ab, bc, ca)
            rlines.add(rline)
            self.play(Write(rline), run_time = 5 / 50)
        self.wait()'''
        lefteqn = TexMobject(
            "\\mathbb{E}(L)"
        ).scale(0.625)
        lefteqn.move_to(4 * UP + 6.5 * LEFT + lefteqn.get_center() - lefteqn.get_critical_point(UL))
        pxeqn = VGroup(
            TexMobject("=\\displaystyle \\frac{\\angle A}{\\pi}\\mathbb{E}(L|S_6,A)+\\frac{\\angle B}{\\pi}\\mathbb{E}(L|S_6,B)+\\frac{\\angle C}{\\pi}\\mathbb{E}(L|S_6,C)").scale(0.625),
            TexMobject("=\\displaystyle \\sum_{x \\in \\{a,b,c\\}}\\frac{h_x}{3\\pi}\\text{sinh}^{-1}\\left(\\frac{x}{p}\\frac{p-x}{s-x}\\right)").scale(0.625),
            #TexMobject("=\\displaystyle \\frac{h_a}{3\\angle A} \\text{sinh}^{-1}\\left(\\frac{a}{p}\\frac{p-a}{s-a}\\right)").scale(0.625),
        )
        pxeqn.arrange(DOWN)
        for obj in pxeqn[1:]:
            obj.align_to(pxeqn[0], LEFT)
        pxeqn.next_to(lefteqn, direction = RIGHT, aligned_edge = UL, buff = 0.9375)
        self.play(Write(lefteqn))
        self.wait()
        for obj in pxeqn:
            self.play(Write(obj))
            self.wait()
        self.play(
            *[FadeOutAndShiftDown(obj) for obj in [lefteqn, pxeqn]]
        )
        self.wait()
        self.wait(5)

def triangle_grp(obj_a, obj_b, obj_c, text_scale = 0.75):
    triangle_color = YELLOW
    altitude_color = GREEN
    a, b, c = obj_a, obj_b, obj_c
    adot, bdot, cdot = Dot(a), Dot(b), Dot(c)
    vertices_grp = VGroup(adot, bdot, cdot)
    alabel, blabel, clabel = TexMobject("A").scale(text_scale), TexMobject("B").scale(text_scale), TexMobject("C").scale(text_scale)
    alabel.next_to(adot, direction = LEFT, buff = 0.1)
    blabel.next_to(bdot, direction = DR, buff = 0.05)
    clabel.next_to(cdot, direction = DL, buff = 0.05)
    ab, bc, ca = Line(a, b), Line(b, c), Line(c, a)
    side_grp = VGroup(ab, bc, ca)
    ddot, edot, fdot = PerpendicularFoot(adot, bc), PerpendicularFoot(bdot, ca), PerpendicularFoot(cdot, ab)
    dlabel, elabel, flabel = TexMobject("D").scale(text_scale), TexMobject("E").scale(text_scale), TexMobject("F").scale(text_scale)
    dlabel.next_to(ddot, direction = DOWN, buff = 0.05)
    elabel.next_to(edot, direction = RIGHT, buff = 0.05)
    flabel.next_to(fdot, direction = LEFT, buff = 0.05)
    label_grp = VGroup(alabel, blabel, clabel, dlabel, elabel, flabel)
    dangle, eangle, fangle = Angle(a, ddot.get_center(), c), Angle(a, edot.get_center(), b), Angle(c, fdot.get_center(), b)
    dangle.scale_in_place(0.5, about_point = ddot.get_center())
    eangle.scale_in_place(0.5, about_point = edot.get_center())
    fangle.scale_in_place(0.5, about_point = fdot.get_center())
    angle_grp = VGroup(dangle, eangle, fangle)
    ad, be, cf = Line(adot.get_center(), ddot.get_center()), Line(bdot.get_center(), edot.get_center()), Line(cdot.get_center(), fdot.get_center())
    altitude_grp = VGroup(ad, be, cf)
    altitude_grp.set_color(altitude_color)
    tri = VMobject(stroke_color = triangle_color)
    tri.set_points_as_corners([a, b, c, a])
    return vertices_grp, side_grp, label_grp, angle_grp, altitude_grp, tri

def triintersector(r_line, obj_ab, obj_bc, obj_ca):
    for k in [obj_ab, obj_bc, obj_ca]:
        tempobj = LineSegIntersection(r_line, k)
        if len(tempobj) > 0:
            return tempobj[0]

def rptangle(obj_a, obj_b, obj_c, randomangle = [0, 2 * PI]):
    ac_rand = random.random()
    ab_rand = random.random()
    if ac_rand + ab_rand > 1:
        ac_rand, ab_rand = 1 - ac_rand, 1 - ab_rand
    rpt = obj_a + (obj_c - obj_a) * ac_rand + (obj_b - obj_a) * ab_rand
    rangle = randomangle[0] + random.random() * (randomangle[1] - randomangle[0])
    return rpt, rangle

def line_creator(pnt, angl, obj_ab, obj_bc, obj_ca):
    rlin = Line(pnt, pnt + 5 * np.array([np.cos(angl), np.sin(angl), 0]))
    pt = triintersector(rlin, obj_ab, obj_bc, obj_ca)
    if pt == None:
        return Line(ORIGIN, ORIGIN + RIGHT / 1000)
    rline = Line(pnt, pt.get_center())
    return rline

def lines_grp(obj_a, obj_b, obj_c, randomangle = [0, 2 * PI], nlimit = 100):
    resgrp = VGroup()
    color_grp = [BLUE, GREEN, TEAL, GOLD, RED, PURPLE, ORANGE, PINK, MAROON]
    leng = len(color_grp)
    cnt = 0
    for _ in range(nlimit):
        tempm, tempn = rptangle(obj_a, obj_b, obj_c, randomangle = randomangle)
        randomline = line_creator(tempm, tempn, Line(obj_a, obj_b), Line(obj_b, obj_c), Line(obj_c, obj_a)).set_color(color_grp[cnt % leng])
        cnt += 1
        resgrp.add(randomline)
    return resgrp

class ProblemIntroduction(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        self.play(
            #self.camera_frame.scale, 1.25,
            self.camera_frame.shift, 2 * UR + RIGHT,
        )
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.play(ShowCreation(triangle), Write(labels[:3]))
        self.wait()
        rlines = VGroup()
        for _ in range(8):
            m, n = rptangle(a, b, c)
        rline = line_creator(m, n, *sides)
        rlines.add(rline)
        self.play(ShowCreation(rline))
        self.wait()
        for _ in range(100):
            m, n = rptangle(a, b, c)
            rline = line_creator(m, n, *sides)
            rlines.add(rline)
            self.play(ShowCreation(rline), run_time = 5 / 100)
        self.wait()
        prereqs = TexMobject("\\text{\\underline{Prerequisite}}").set_color(YELLOW)
        prereqs.shift(2 * UR + RIGHT + 3.25 * UP)
        self.play(FadeInFrom(prereqs), FadeOut(triangle), FadeOut(rlines), FadeOut(labels[:3]))
        self.wait()
        reqpts = VGroup(
            TextMobject("""$\\bullet$ Given a random variable $X$ with pdf $f(x)$, $\\mathbb{E}(g(X))=\\displaystyle\\int\\limits_{-\\infty}^\\infty g(x)f(x)\\,dx$"""),
            TextMobject("""$\\circ$ If $X \\sim U(0, 2\\pi)$, then $\\mathbb{E}(g(X))=\\displaystyle\\int\\limits_{0}^{2\\pi} g(x)\\,\\frac{dx}{2\\pi}$""", tex_to_color_map = {"$X \\sim U(0, 2\\pi)$": BLUE}),
            TextMobject("""$\\circ$ If $P(x,y) \\sim \\triangle ABC$, then $\\mathbb{E}(g(P))=\\displaystyle\\iint\\limits_{(x,y)\\in\\triangle ABC} g(P)\\,\\frac{dydx}{\\triangle ABC}$""", tex_to_color_map = {"$P(x,y) \\sim \\triangle ABC$": BLUE}),
            TextMobject("""$\\bullet$ $\\displaystyle\\int\\limits_{0}^{\\tan^{-1}(x)} \\text{sec}(\\theta) \\,d\\theta=\\sinh^{-1}(x)$"""),
            TextMobject("""$\\bullet$ $\\mathbb{E}(Y)=\\mathbb{E}(\\mathbb{E}(Y|X))$ (or) $\\mathbb{E}(Y)=\\displaystyle\\int\\limits_{-\\infty}^\\infty \\mathbb{E}(Y|X=x)f(x)\\,dx$"""),
        )
        for text in reqpts:
            text.scale(0.75)
        reqpts.arrange(DOWN, buff = 0.05)
        for text in reqpts[1:]:
            text.align_to(reqpts[0], LEFT)
        reqpts.move_to(3 * LEFT + 4.75 * UP + reqpts.get_center() - reqpts.get_critical_point(UL))
        for text in reqpts:
            self.play(Write(text))
            self.wait()
        self.play(
            *[FadeOut(obj) for obj in [prereqs, reqpts]],
            *[FadeIn(obj) for obj in [triangle, labels[:3]]]
        )
        self.wait()
        self.play(Write(coordsys))
        self.wait()
        self.play(*[Write(obj) for obj in [labels[3:], angles, altitudes]])
        self.wait()
        altitude_lengths = VGroup(
            TexMobject("|AD|=h_a"),
            TexMobject("|BE|=h_b"),
            TexMobject("|CF|=h_c"),
        )
        for obj in altitude_lengths:
            obj.scale_in_place(0.75)
        altitude_lengths.arrange(DOWN, center = False)
        altitude_lengths.move_to(6 * RIGHT + 3 * UP)
        self.play(Write(altitude_lengths))
        self.wait()
        rline = line_creator(UR, 35 * DEGREES, *sides)
        self.play(
            *[FadeOut(obj) for obj in [labels[3:], angles, altitudes, altitude_lengths]],
            self.camera_frame.shift, 2 * RIGHT,
            ShowCreation(rline)
        )
        self.wait()
        xylabel = TexMobject("P(x,y)").scale(0.75).next_to(rline.get_start(), DOWN, buff = 0.05)
        lengthlabel = TexMobject("L").scale(0.75).next_to(rline.get_center(), UP, buff = 0.2)
        angledashedline = DashedLine(rline.get_start(), rline.get_start() + RIGHT).set_color(GREY)
        angle = Angle(angledashedline.get_end(), rline.get_start(), rline.get_end())
        tlabel = TexMobject("t").scale(0.75).next_to(angledashedline.get_center(), RIGHT, buff = 0.1).shift(0.225 * UP)
        self.play(*[Write(obj) for obj in (xylabel, lengthlabel, angledashedline, angle, tlabel)])
        eqnlist = VGroup(
            TextMobject("$L$ is some function $g$ of $x$, $y$ and $t$"),
            TexMobject("\\mathbb{E}(L)=\\displaystyle \\int\\limits_0^{2\\pi}\\int g(x,y,t)\\,\\frac{dydx}{\\triangle ABC}\\frac{dt}{2\\pi}", tex_to_color_map = {"g(x,y,t)": RED}),
            TextMobject("But what is $g(x,y,t)$?", tex_to_color_map = {"$g(x,y,t)$": RED}),
        )
        for text in eqnlist:
            text.scale(0.75)
        eqnlist.arrange(DOWN, center = True)
        eqnlist.move_to(5 * RIGHT + 4 * UP + eqnlist.get_center() - eqnlist.get_critical_point(UL))
        self.play(Write(eqnlist[:2]))
        self.wait()
        self.play(Write(eqnlist[-1]))
        self.wait()
        cansolvetxt = TextMobject("Can this be solved \\\\ in closed form?")
        tbubble = ThoughtBubble()
        tbubble.pin_to(eqnlist[2])
        tbubble.add_content(cansolvetxt)
        tbubble.resize_to_content()
        bub_grp = VGroup(tbubble, cansolvetxt)
        self.play(Write(bub_grp))
        self.wait()
        self.play(
            *[FadeOutAndShiftDown(obj) for obj in [rline, eqnlist, bub_grp, xylabel, lengthlabel, angledashedline, angle, tlabel]],
            self.camera_frame.shift, UR,
        )
        '''testeqn = TexMobject(
            "\\mathbb{E}(L|S_6,A)=\\displaystyle \\frac{\\tan^{-1}(CD/AD)}{\\angle A}\\mathbb{E}(L|S_5,CD)+\\frac{\\tan^{-1}(BD/AD)}{\\angle A}\\mathbb{E}(L|S_5,BD)"
        ).scale(0.625)
        testeqn.move_to(8 * RIGHT + 4 * UP)
        self.play(Write(testeqn))'''
        self.wait(5)

class ProblemSetup1(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        self.add(coordsys, triangle, labels[:3])
        conditioningtext = TextMobject("Conditioning is the soul of Statistics", tex_to_color_map = {"Conditioning": YELLOW})
        conditioningtext.shift(3 * UR + 3 * RIGHT + 3 * UP)
        self.play(Write(conditioningtext))
        self.wait()
        self.play(FadeOutAndShift(conditioningtext, DOWN))
        self.wait()
        setupdescription = TextMobject("Setup S1: Choose a point randomly in $\\triangle ABC$ and draw a line perpendicular to base $BC$ from that point.", tex_to_color_map = {"Setup S1": YELLOW}).scale(0.75)
        setupdescription.shift(3 * UR + 3 * RIGHT + 3.25 * UP)
        self.play(Write(setupdescription))
        rlines = lines_grp(a, b, c, randomangle = [-PI / 2 - 1 / 1000, -PI / 2 + 1 / 1000])
        for obj in rlines[:50]:
            self.play(ShowCreation(obj), run_time = 5 / 50)
        self.wait()
        self.play(FadeOut(rlines[:50]), ShowCreation(rlines[79]))
        self.wait()
        xylabel = TexMobject("P(x,y)").scale(0.75).next_to(rlines[79].get_start(), RIGHT, buff = 0.125)
        self.play(Write(xylabel))
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}(L|S_1)").scale(0.75)
        equationslist = VGroup(
            TexMobject("=\\displaystyle \\int\\limits_{(x,y)\\in \\triangle ABC}y\\,\\frac{dydx}{\\triangle ABC}"),
            TexMobject("=\\displaystyle \\frac{|AD|}{3}")
        )
        for text in equationslist:
            text.scale(0.75)
        equationslist.arrange(DOWN)
        for text in equationslist[1:]:
            text.align_to(equationslist[0], LEFT)
        expectedeqn.move_to(6 * RIGHT + 4 * UP)
        equationslist.move_to(expectedeqn.get_critical_point(RIGHT) + equationslist.get_center() - equationslist[0][0][0].get_critical_point(LEFT) + 0.25 * RIGHT)
        #equationslist.next_to(expectedeqn)
        #self.play(Write(equationslist[0][0][0].copy().scale(5)))
        #print(expectedeqn.get_center())
        self.play(Write(expectedeqn), Write(equationslist[0]))
        self.wait()
        self.play(
            ShowCreation(altitudes[0]),
            Write(labels[3]),
            Write(equationslist[1])
        )
        self.wait()
        self.play(*[FadeOut(obj) for obj in [setupdescription, expectedeqn, equationslist, rlines[79], xylabel]])
        self.wait(5)

class ProblemSetup2(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        self.add(coordsys, triangle, labels[:3], labels[3], altitudes[0])
        self.wait()
        setupdescription = TextMobject("Setup S2: Given $M$ between $D$ and $C$, choose a point randomly in $\\triangle ABC$ and draw a line parallel to $AM$ $towards$ base $BC$.", tex_to_color_map = {"Setup S2": YELLOW}).scale(0.75)
        setupdescription.shift(3 * UR + 3 * RIGHT + 3.25 * UP)
        mdot = Dot(3 * RIGHT)
        mlabel = TexMobject("M").scale(0.75)
        mlabel.next_to(mdot, DOWN, 0.05)
        am = DashedLine(a, mdot.get_center(), dash_length = 0.125).set_color(GREY)
        self.play(
            Write(setupdescription),
            *[Write(obj) for obj in [mlabel, am]]
        )
        self.wait()
        fixedangle = am.get_angle()
        rlines = lines_grp(a, b, c, randomangle = [fixedangle - 1 / 1000, fixedangle + 1 / 1000])
        for obj in rlines[:50]:
            self.play(ShowCreation(obj), run_time = 5 / 50)
        self.wait()
        self.play(FadeOut(rlines[:50]), ShowCreation(rlines[79]))
        self.wait()
        xylabel = TexMobject("P(x,y)").scale(0.75).next_to(rlines[79].get_start(), LEFT, buff = 0.125)
        qlabel = TexMobject("Q").scale(0.75).next_to(rlines[79].get_end(), DOWN, buff = 0.125)
        self.play(Write(xylabel), Write(qlabel))
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}_{ABC}(L|S_2,AM)").scale(0.625)
        equationslist = VGroup(
            TexMobject("=\\displaystyle \\int\\limits_{(x,y)\\in \\triangle ABC} |PQ| \\,\\frac{dydx}{\\triangle ABC}"),
            TexMobject("=\\displaystyle \\int\\limits_{(x,y)\\in \\triangle ABC} \\frac{y}{\\cos t} \\,\\frac{dydx}{\\triangle ABC}"),
            TexMobject("=\\displaystyle \\frac{|AD|}{3\\cos t}=\\frac{|AM|}{3}")
        )
        for text in equationslist:
            text.scale(0.625)
        equationslist.arrange(DOWN)
        for text in equationslist[1:]:
            text.align_to(equationslist[0], LEFT)
        expectedeqn.move_to(6 * RIGHT + 4 * UP)
        equationslist.move_to(expectedeqn.get_critical_point(RIGHT) + equationslist.get_center() - equationslist[0][0][0].get_critical_point(LEFT) + 0.25 * RIGHT)
        self.play(Write(expectedeqn), Write(equationslist[0]))
        self.wait()
        tangle = Angle(mdot.get_center(), a, altitudes[0].get_end(), radius = 0.75)
        tlabel = TexMobject("t").scale(0.75).move_to(1.75 * RIGHT + 3 * UP + 0.0625 * LEFT)
        tgroup = VGroup(tangle, tlabel)
        self.play(Write(tgroup))
        self.wait()
        for text in equationslist[1:]:
            self.play(Write(text))
            self.wait()
        self.play(*[FadeOut(obj) for obj in [setupdescription, expectedeqn, equationslist, rlines[79], xylabel, qlabel]])
        self.wait(5)

class ProblemSetup3(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        setupdescription = TextMobject("Setup S3: Given $M$ between $D$ and $C$, choose a point randomly in $\\triangle ABC$ and draw a line parallel to $AM$ $away$ from base $BC$.", tex_to_color_map = {"Setup S3": YELLOW}).scale(0.75)
        setupdescription.shift(3 * UR + 3 * RIGHT + 3.25 * UP)
        mdot = Dot(3 * RIGHT)
        mlabel = TexMobject("M").scale(0.75)
        mlabel.next_to(mdot, DOWN, 0.05)
        am = DashedLine(a, mdot.get_center(), dash_length = 0.125).set_color(GREY)
        tangle = Angle(mdot.get_center(), a, altitudes[0].get_end(), radius = 0.75)
        tlabel = TexMobject("t").scale(0.75).move_to(1.75 * RIGHT + 3 * UP + 0.0625 * LEFT)
        tgroup = VGroup(tangle, tlabel)
        self.add(coordsys, triangle, labels[:3], mlabel, labels[3], altitudes[0], am, tgroup)
        self.wait()
        self.play(
            Write(setupdescription),
        )
        self.wait()
        fixedangle = PI + am.get_angle()
        rlines = lines_grp(a, b, c, randomangle = [fixedangle - 1 / 1000, fixedangle + 1 / 1000])
        for obj in rlines[:50]:
            self.play(ShowCreation(obj), run_time = 5 / 50)
        self.wait()
        lineframe = 63
        self.play(FadeOut(rlines[:50]), ShowCreation(rlines[lineframe]))
        self.wait()
        xylabel = TexMobject("P(x,y)").scale(0.75).next_to(rlines[lineframe].get_start(), RIGHT, buff = 0.125)
        self.play(
            Write(xylabel),
        )
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}_{ABC}(L|S_3,AM)").scale(0.5625)
        equationslist = VGroup(
            TexMobject("=\\displaystyle \\frac{\\triangle MAB}{\\triangle ABC}\\mathbb{E}(L|P \\in \\triangle MAB) + \\frac{\\triangle MCA}{\\triangle ABC}\\mathbb{E}(L|P \\in \\triangle MCA)"),
            TexMobject("=\\displaystyle \\frac{\\triangle MAB}{\\triangle ABC}\\mathbb{E}_{MAB}(L|S_2,MA) + \\frac{\\triangle MCA}{\\triangle ABC}\\mathbb{E}_{MCA}(L|S_2,MA)"),
            TexMobject("=\\displaystyle \\frac{\\triangle MAB}{\\triangle ABC}\\frac{|MA|}{3} + \\frac{\\triangle MCA}{\\triangle ABC}\\frac{|MA|}{3}"),
            TexMobject("=\\displaystyle \\frac{|AM|}{3}"),
        )
        for text in equationslist:
            text.scale(0.5625)
        equationslist.arrange(DOWN)
        for text in equationslist[1:]:
            text.align_to(equationslist[0], LEFT)
        expectedeqn.move_to(4 * RIGHT + 4 * UP)
        equationslist.move_to(expectedeqn.get_critical_point(RIGHT) + equationslist.get_center() - equationslist[0][0][0].get_critical_point(LEFT) + 0.25 * RIGHT)
        self.play(Write(expectedeqn), Write(equationslist[0]))
        self.wait()
        linesfrommca = lines_grp(mdot.get_center(), c, a, randomangle = [fixedangle - 1 / 1000, fixedangle + 1 / 1000])
        for obj in linesfrommca[:50]:
            self.play(Write(obj), run_time = 5 / 100)
        self.wait()
        self.play(FadeOut(linesfrommca[:50]))
        self.wait()
        for text in equationslist[1:]:
            self.play(Write(text))
            self.wait()
        self.play(*[FadeOut(obj) for obj in [setupdescription, expectedeqn, equationslist, rlines[lineframe], xylabel]])
        #self.play(*[FadeOut(obj) for obj in [tgroup, mlabel, am, altitudes[0], labels[3]]])
        self.wait(5)

class ProblemSetup4(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        setupdescription = TextMobject("Setup S4: Given $M$ between $D$ and $C$, choose a point randomly in $\\triangle ABC$ and draw a line parallel to $AM$.", tex_to_color_map = {"Setup S4": YELLOW}).scale(0.75)
        setupdescription.shift(3 * UR + 3 * RIGHT + 3.25 * UP)
        mdot = Dot(3 * RIGHT)
        mlabel = TexMobject("M").scale(0.75)
        mlabel.next_to(mdot, DOWN, 0.05)
        am = DashedLine(a, mdot.get_center(), dash_length = 0.125).set_color(GREY)
        tangle = Angle(mdot.get_center(), a, altitudes[0].get_end(), radius = 0.75)
        tlabel = TexMobject("t").scale(0.75).move_to(1.75 * RIGHT + 3 * UP + 0.0625 * LEFT)
        tgroup = VGroup(tangle, tlabel)
        self.add(coordsys, triangle, labels[:3], mlabel, labels[3], altitudes[0], am, tgroup)
        self.wait()
        self.play(
            Write(setupdescription),
        )
        self.wait()
        fixedangle1 = PI + am.get_angle()
        fixedangle2 = am.get_angle()
        rlines1 = lines_grp(a, b, c, randomangle = [fixedangle1 - 1 / 1000, fixedangle1 + 1 / 1000])
        rlines2 = lines_grp(a, b, c, randomangle = [fixedangle2 - 1 / 1000, fixedangle2 + 1 / 1000])
        rlines = VGroup()
        for i in range(100):
            if i & 1:
                rlines.add(rlines1[i])
            else:
                rlines.add(rlines2[i])
        for obj in rlines[:50]:
            self.play(ShowCreation(obj), run_time = 5 / 50)
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}(L|S_4,AM)").scale(0.625)
        equationslist = VGroup(
            TextMobject("$=p\\text{ }\\mathbb{E}(L|$towards BC$) + (1-p)\\text{ }\\mathbb{E}(L|$away from BC$)$", tex_to_color_map = {"towards BC": BLUE, "away from BC": PURPLE}),
            TextMobject("$=\\displaystyle p\\text{ }\\mathbb{E}_{ABC}(L|$$S_2$$,AM) + (1-p)\\text{ }\\mathbb{E}_{ABC}(L|$$S_3$$,AM)$", tex_to_color_map = {"$S_2$": BLUE, "$S_3$": PURPLE}),
            TexMobject("=\\displaystyle p\\text{ }\\frac{|AM|}{3} + (1-p)\\text{ }\\frac{|AM|}{3}"),
            TexMobject("=\\displaystyle \\frac{|AM|}{3}"),
        )
        for text in equationslist:
            text.scale(0.625)
        equationslist.arrange(DOWN)
        for text in equationslist[1:]:
            text.align_to(equationslist[0], LEFT)
        expectedeqn.move_to(4.5 * RIGHT + 4 * UP)
        equationslist.move_to(expectedeqn.get_critical_point(RIGHT) + equationslist.get_center() - equationslist[0][0][0].get_critical_point(LEFT) + 0.125 * RIGHT)
        self.play(Write(expectedeqn), Write(equationslist[0]))
        self.wait()
        for text in equationslist[1:]:
            self.play(Write(text))
            self.wait()
        self.play(*[FadeOut(obj) for obj in [expectedeqn, equationslist, rlines[:50]]])
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}(L|S_4,t)=\\displaystyle\\frac{|AM|}{3}=\\frac{|AD|}{3\\cos t}=\\frac{h_a}{3\\cos t}").scale(0.75)
        expectedeqn.move_to(4.5 * RIGHT + 4 * UP + expectedeqn.get_center() - expectedeqn.get_critical_point(UL))
        self.play(Write(expectedeqn))
        self.wait()
        self.play(*[FadeOut(obj) for obj in [setupdescription, expectedeqn, tgroup, am, mlabel]])
        #self.play(*[FadeOut(obj) for obj in [tgroup, mlabel, am, altitudes[0], labels[3]]])
        self.wait(5)

class ProblemSetup5(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        setupdescription = TextMobject("Setup S5: Choose a point randomly in $\\triangle ABC$ and draw a line (towards or away from $BC$) such that the angle $t$ from $AD$ is $U(0,\\tan^{-1}(CD/AD))$.", tex_to_color_map = {"Setup S5": YELLOW}).scale(0.75)
        setupdescription.shift(3 * UR + 3 * RIGHT + 3.25 * UP)
        mval = ValueTracker(3)
        mdot = Dot(mval.get_value() * RIGHT)
        mlabel = TexMobject("M").scale(0.75)
        mlabel.next_to(mdot, DOWN, 0.05)
        am = DashedLine(a, mdot.get_center(), dash_length = 0.125).set_color(GREY)
        tangle = Angle(altitudes[0].get_end(), a, c, radius = 0.75)
        tlabel = TexMobject("t").scale(0.75).move_to(1.75 * RIGHT + 3 * UP + 0.0625 * LEFT)
        tgroup = VGroup(tangle, tlabel)
        self.add(coordsys, triangle, labels[:3], labels[3], altitudes[0])
        self.wait()
        self.play(
            Write(setupdescription),
            ShowCreation(tangle),
        )
        rlines1 = lines_grp(a, b, c, randomangle = [-PI / 2 , am.get_angle()])
        rlines2 = lines_grp(a, b, c, randomangle = [PI / 2, PI + am.get_angle()])
        rlines = VGroup()
        for i in range(100):
            if i & 1:
                rlines.add(rlines1[i])
            else:
                rlines.add(rlines2[i])
        for obj in rlines[:50]:
            self.play(ShowCreation(obj), run_time = 5 / 50)
        self.wait()
        self.play(FadeOut(rlines[:50]))
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}(L|S_5,CD)").scale(0.625)
        equationslist = VGroup(
            TexMobject("=\\displaystyle \\mathbb{E}(\\mathbb{E}(L|S_4,t))"),
            TexMobject("=\\displaystyle \\frac{1}{\\tan^{-1}(CD/AD)}\\int\\limits_{0}^{\\tan^{-1}(CD/AD)}\\frac{h_a}{3\\cos t}\\,dt"),
            TexMobject("=\\displaystyle \\frac{h_a}{3}\\frac{\\text{sinh}^{-1}(CD/AD)}{\\tan^{-1}(CD/AD)}"),
        )
        for text in equationslist:
            text.scale(0.625)
        equationslist.arrange(DOWN)
        for text in equationslist[1:]:
            text.align_to(equationslist[0], LEFT)
        expectedeqn.move_to(5 * RIGHT + 4 * UP)
        equationslist.move_to(expectedeqn.get_critical_point(RIGHT) + equationslist.get_center() - equationslist[0][0][0].get_critical_point(LEFT) + 0.125 * RIGHT)
        self.play(Write(expectedeqn), Write(equationslist[0]))
        self.wait()
        for text in equationslist[1:]:
            self.play(Write(text))
            self.wait()
        self.play(*[FadeOut(obj) for obj in [tangle, setupdescription, expectedeqn, equationslist]])
        #self.play(*[FadeOut(obj) for obj in [tgroup, mlabel, am, altitudes[0], labels[3]]])
        self.wait(5)

class ProblemSetup6(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        setupdescription = TextMobject("Setup S6: Choose a point randomly in $\\triangle ABC$ and draw a line (towards or away from $BC$) such that it's angle $t$ from $AC$ is $U(0,\\angle A)$.", tex_to_color_map = {"Setup S6": YELLOW}).scale(0.75)
        setupdescription.shift(3 * UR + 3 * RIGHT + 3.25 * UP)
        tangle = Angle(c, a, b, radius = 0.75)
        tlabel = TexMobject("t").scale(0.75).move_to(1.75 * RIGHT + 3 * UP + 0.0625 * LEFT)
        tgroup = VGroup(tangle, tlabel)
        self.add(coordsys, triangle, labels[:3], labels[3], altitudes[0])
        self.wait()
        self.play(
            Write(setupdescription),
            ShowCreation(tangle),
        )
        rlines1 = lines_grp(a, b, c, randomangle = [2 * PI + sides[0].get_angle() - PI, sides[2].get_angle()])
        rlines2 = lines_grp(a, b, c, randomangle = [2 * PI + sides[0].get_angle(), PI + sides[2].get_angle()])
        rlines = VGroup()
        for i in range(100):
            if i & 1:
                rlines.add(rlines1[i])
            else:
                rlines.add(rlines2[i])
        for obj in rlines[:50]:
            self.play(ShowCreation(obj), run_time = 5 / 50)
        self.wait()
        self.play(FadeOut(rlines[:50]))
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}(L|S_6,A)").scale(0.625)
        equationslist = VGroup(
            TexMobject("=\\displaystyle \\frac{\\angle DAB}{\\angle A}\\mathbb{E}(L|S_5,BD)+\\frac{\\angle CAD}{\\angle A}\\mathbb{E}(L|S_5,CD)"),
            TexMobject("=\\displaystyle \\frac{h_a}{3\\angle A}\\text{sinh}^{-1}\\left(\\frac{BD}{AD}\\right)+\\frac{h_a}{3\\angle A}\\text{sinh}^{-1}\\left(\\frac{CD}{AD}\\right)"),
            TexMobject("=\\displaystyle \\frac{h_a}{3\\angle A}\\text{sinh}^{-1}\\left(\\frac{a}{p}\\frac{p-a}{s-a}\\right)"),
            TextMobject("where $p=2s=a+b+c$"),
        )
        for text in equationslist:
            text.scale(0.625)
        equationslist.arrange(DOWN)
        for text in equationslist[1:]:
            text.align_to(equationslist[0], LEFT)
        expectedeqn.move_to(5 * RIGHT + 4 * UP)
        equationslist.move_to(expectedeqn.get_critical_point(RIGHT) + equationslist.get_center() - equationslist[0][0][0].get_critical_point(LEFT) + 0.125 * RIGHT)
        self.play(Write(expectedeqn), Write(equationslist[0]))
        self.wait()
        leftangle = Angle(b, a, altitudes[0].get_end(), radius = 0.75, color = BLUE)
        rightangle = Angle(c, a, altitudes[0].get_end(), radius = 0.75, color = PURPLE)
        self.play(
            FadeOut(tangle),
            ShowCreation(leftangle),
            ShowCreation(rightangle)
        )
        self.wait()
        for text in equationslist[1:]:
            self.play(Write(text))
            self.wait()
        self.play(*[FadeOut(obj) for obj in [labels[3], altitudes[0], leftangle, rightangle, setupdescription, expectedeqn, equationslist]])
        #self.play(*[FadeOut(obj) for obj in [tgroup, mlabel, am, altitudes[0], labels[3]]])
        self.wait(5)

class ProblemSetup7(GraphScene, MovingCameraScene):
    def setup(self):
        GraphScene.setup(self)
        MovingCameraScene.setup(self)
    
    def construct(self):
        self.setup_axes()
        self.play(FadeOut(self.axes))
        self.wait()
        nplane = VGroup(DoubleArrow(1 * LEFT, 7 * RIGHT), DoubleArrow(1 * DOWN, 4 * UP))
        coord_labels = VGroup(TexMobject("x"), TexMobject("y"))
        coord_labels[0].next_to(nplane[0].get_end(), direction = DOWN, buff = 0.25)
        coord_labels[1].next_to(nplane[1].get_end(), direction = RIGHT, buff = 0.25)
        coordsys = VGroup(nplane, coord_labels)
        a, b, c = 4 * UP + 1.5 * RIGHT, ORIGIN, 6 * RIGHT
        verts, sides, labels, angles, altitudes, triangle = triangle_grp(a, b, c)
        self.camera_frame.shift(3 * UR + 3 * RIGHT)
        setupdescription = TextMobject("Choose a point randomly in $\\triangle ABC$ and draw a line at a random angle from that point until it meets the triangle.").scale(0.75)
        setupdescription.shift(3 * UR + 3 * RIGHT + 3.25 * UP)
        tangle = Angle(c, a, b, radius = 0.75)
        tlabel = TexMobject("t").scale(0.75).move_to(1.75 * RIGHT + 3 * UP + 0.0625 * LEFT)
        tgroup = VGroup(tangle, tlabel)
        self.add(coordsys, triangle, labels[:3])
        self.wait()
        self.play(
            Write(setupdescription),
        )
        rlines = lines_grp(a, b, c)
        for obj in rlines[:50]:
            self.play(ShowCreation(obj), run_time = 5 / 50)
        self.wait()
        self.play(FadeOut(rlines[:50]))
        self.wait()
        expectedeqn = TexMobject("\\mathbb{E}(L)").scale(0.625)
        equationslist = VGroup(
            TexMobject("=\\displaystyle \\frac{\\angle A}{\\pi}\\mathbb{E}(L|S_6,A)+\\frac{\\angle B}{\\pi}\\mathbb{E}(L|S_6,B)+\\frac{\\angle C}{\\pi}\\mathbb{E}(L|S_6,C)"),
            TexMobject("=\\displaystyle \\sum_{k \\in \\{a,b,c\\}}\\frac{h_k}{3\\pi}\\text{sinh}^{-1}\\left(\\frac{k}{p}\\frac{p-k}{s-k}\\right)"),
        )
        for text in equationslist:
            text.scale(0.625)
        equationslist.arrange(DOWN)
        for text in equationslist[1:]:
            text.align_to(equationslist[0], LEFT)
        expectedeqn.move_to(5 * RIGHT + 4 * UP)
        equationslist.move_to(expectedeqn.get_critical_point(RIGHT) + equationslist.get_center() - equationslist[0][0][0].get_critical_point(LEFT) + 0.125 * RIGHT)
        self.play(Write(expectedeqn), Write(equationslist[0]))
        self.wait()
        leftangle = Angle(a, b, c, radius = 0.75, color = BLUE)
        rightangle = Angle(b, c, a, radius = 0.75, color = PURPLE)
        self.play(
            ShowCreation(tangle),
            ShowCreation(leftangle),
            ShowCreation(rightangle)
        )
        for text in equationslist[1:]:
            self.play(Write(text))
            self.wait()
        finaleqn = TexMobject("\\displaystyle \\mathbb{E}(L)=\\frac{h_a}{3\\pi}\\text{sinh}^{-1}\\left(\\frac{a}{p}\\frac{p-a}{s-a}\\right)+\\frac{h_b}{3\\pi}\\text{sinh}^{-1}\\left(\\frac{b}{p}\\frac{p-b}{s-b}\\right)+\\frac{h_c}{3\\pi}\\text{sinh}^{-1}\\left(\\frac{c}{p}\\frac{p-c}{s-c}\\right)").scale(0.75)
        finaleqn.move_to(triangle.get_center() + 3.5 * DOWN)
        self.play(
            *[FadeOut(obj) for obj in [setupdescription, expectedeqn, equationslist, coordsys]],
            self.camera_frame.move_to, 2 * UR + DOWN + RIGHT,
            FadeIn(finaleqn),
        )
        self.wait()
        #self.play(*[FadeOut(obj) for obj in [tgroup, mlabel, am, altitudes[0], labels[3]]])
        self.wait(5)

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    #command_B = module_name + " " + "ProblemSetup7" + " -pl"
    #command_B = module_name + " " + "IntroductionScene" + " -p",
    command_B = module_name + " -p"
    os.system(clear_cmd)
    os.system(command_A + command_B)