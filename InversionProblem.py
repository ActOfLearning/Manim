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

class RoseCurve(ParametricFunction):
    CONFIG = {
        "radius"        : 1,
        "start_theta"   : 0,
        "kval"          : 3 / 2,
        "end_theta"     : 2 * PI,
        "step_size"     : 0.001,
        "density"       : 50 * DEFAULT_POINT_DENSITY_1D,
        "color"         : BLUE,
    }
    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        ParametricFunction.__init__(self, self.pos_func, **kwargs)

    def pos_func(self, t):
        T = self.start_theta + t * (self.end_theta - self.start_theta)
        return self.radius * np.array([
            np.cos(self.kval * T) * np.cos(T),
            np.cos(self.kval * T) * np.sin(T),
            0
        ])

class ProblemIntroduction(GraphScene, MovingCameraScene):
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
        def Circle_Inverse_Func(grap, x_min = -10, x_max = 10, step_size = 0.01, about_point = ORIGIN, inversion_radius = 1):
            res_vm = VMobject()
            res_list = [Circle_Inversion(Dot(self.input_to_graph_point(x, grap)), about_point, inversion_radius).get_center() for x in np.arange(x_min, x_max, step_size)]
            res_vm.set_points_as_corners(res_list)
            return res_vm
        #ra, rb, rc = 1, 1.4, 1.6
        ra, rb, rc = 1, 1.6, 1.9
        ca = Circle(radius = ra).move_to(ra * LEFT).set_color(self.given_circle_color)
        cb = Circle(radius = rb).move_to(rb * RIGHT).set_color(self.given_circle_color)
        xc = rc * (ra - rb) / (ra + rb)
        yc = ra * rb * rc * (ra + rb + rc) / (ra + rb) ** 2
        yc = 2 * yc ** 0.5
        cena = -ra
        cenb = rb
        cenc = xc + yc * 1j
        curva, curvb, curvc = 1 / ra, 1 / rb, 1 / rc
        curvi = curva + curvb + curvc + 2 * (curva * curvb + curvb * curvc + curvc * curva) ** 0.5
        curvo = curva + curvb + curvc - 2 * (curva * curvb + curvb * curvc + curvc * curva) ** 0.5
        ri, ro = abs(1 / curvi), abs(1 / curvo)
        ceni = curva * cena + curvb * cenb + curvc * cenc + 2  * (curva * cena * curvb * cenb + curvb * cenb * curvc * cenc + curvc * cenc * curva * cena) ** 0.5
        ceni /= curvi
        ceno = curva * cena + curvb * cenb + curvc * cenc - 2  * (curva * cena * curvb * cenb + curvb * cenb * curvc * cenc + curvc * cenc * curva * cena) ** 0.5
        ceno /= curvo
        cena *= RIGHT
        cenb *= RIGHT
        cenc = cenc.real * RIGHT + cenc.imag * UP
        ceni = ceni.real * RIGHT + ceni.imag * UP
        ceno = ceno.real * RIGHT + ceno.imag * UP
        ri, ro = abs(1 / curvi), abs(1 / curvo)
        ci = Circle(radius = ri).move_to(ceni).set_color(self.incircle_color)
        co = Circle(radius = ro).move_to(ceno).set_color(self.excircle_color)
        cc = Circle(radius = rc).move_to(cenc).set_color(self.given_circle_color)
        cid, cod = Dot(ceni).scale(0.5).set_color(self.incircle_color), Dot(ceno).scale(0.5).set_color(self.excircle_color)
        setup_grp = VGroup(ca, cb, cc, ci, co, cid, cod)
        setup_grp.shift(-ceno)
        #self.camera_frame.shift(0.25 * UP)
        #titl = TexMobject("\\text{\\underline{Here's a nice puzzle...}}").move_to(ro * UP + 0.75 * UP)
        labela = TextMobject("A").scale(0.5).move_to(ca.get_center())
        labelb = TextMobject("B").scale(0.5).move_to(cb.get_center())
        labelc = TextMobject("C").scale(0.5).move_to(cc.get_center())
        #self.play(Write(titl))
        #self.wait()
        self.play(Write(ca), Write(cb), Write(cc))
        self.play(Write(labela), Write(labelb), Write(labelc))
        self.wait()
        self.play(Write(co), Write(cod))
        self.wait()
        self.play(Write(ci), Write(cid))
        self.wait()
        #dist_line = Line(cid.get_center(), cod.get_center()).set_color_by_gradient([self.incircle_color, self.excircle_color])
        dist_line = Line(cid.get_center(), cod.get_center()).set_color(BLUE)
        fade_val = 0.75
        self.play(
            Write(dist_line), 
            ca.fade, fade_val,
            cb.fade, fade_val,
            cc.fade, fade_val,
            labela.fade, fade_val,
            labelb.fade, fade_val,
            labelc.fade, fade_val,
        )
        self.wait()
        desc_thm = TexMobject("(r_a^{-1}+r_b^{-1}+r_c^{-1}+r^{-1})^2=2(r_a^{-2}+r_b^{-2}+r_c^{-2}+r^{-2})").scale(0.8).move_to(rc * DOWN + 2.25 * DOWN)
        desc_thm_tmp = TexMobject("=r_a^{-1}+r_b^{-1}+r_c^{-1}\\pm 2\\sqrt{r_a^{-1}r_b^{-1}+r_b^{-1}r_c^{-1}+r_c^{-1}r_a^{-1}}").scale(0.8).move_to(rc * DOWN + 2.25 * DOWN + 0.5 * RIGHT)
        rinv = TexMobject("r^{-1}").scale(0.8).next_to(desc_thm_tmp, direction = LEFT, buff = 0.1)
        rlabel = TexMobject("r^{-1},R^{-1}", tex_to_color_map = {"r": self.incircle_color, "R": self.excircle_color, "d": BLUE}).scale(0.8).next_to(desc_thm_tmp, direction = LEFT, buff = 0.1)
        #rlabel[0][0].set_color(self.incircle_color)
        #rlabel[0][3].set_color(self.incircle_color)
        #rlabel[0][4].set_color(self.excircle_color)
        #rlabel[0][8].set_color(self.excircle_color)
        desc_thm_gp = VGroup(rinv, desc_thm_tmp)
        self.play(dist_line.fade, fade_val, Write(desc_thm), self.camera_frame.shift, 0.75 * DOWN)
        self.wait()
        self.play(Transform(desc_thm, desc_thm_gp))
        self.wait()
        self.play(Transform(desc_thm[0], rlabel))
        self.wait()
        reverse_fade = 1 - 1 / (1 - fade_val)
        self.play(
            self.camera_frame.shift, 0.75 * UP,
            ci.fade, fade_val,
            co.fade, fade_val,
            cod.fade, fade_val,
            cid.fade, fade_val,
            #dist_line.fade, fade_val,
            ca.fade, reverse_fade,
            cb.fade, reverse_fade,
            cc.fade, reverse_fade,
            labela.fade, reverse_fade,
            labelb.fade, reverse_fade,
            labelc.fade, reverse_fade,
            FadeOut(desc_thm)
        )
        self.wait()
        self.play(
            ci.fade, reverse_fade,
            co.fade, reverse_fade,
            cod.fade, reverse_fade,
            cid.fade, reverse_fade,
            dist_line.fade, reverse_fade,
        )
        self.wait()
        self.play(WiggleOutThenIn(dist_line, scale_value = 2),)
        self.wait(3)
        #self.remove(titl)
        #titl_temp = TexMobject("\\text{\\underline{Choosing a right inversion circle...}}").move_to(titl.get_center())
        #self.play(Transform(titl, titl_temp))
        #self.wait()
        t = dist_line.get_angle()
        def rot(obj):
            obj.rotate_about_origin(-t)
            obj.rotate_in_place(t)
            return obj
        self.play(
            ci.rotate_about_origin, -t,
            ca.rotate_about_origin, -t,
            cb.rotate_about_origin, -t,
            cc.rotate_about_origin, -t,
            cid.rotate_about_origin, -t,
            dist_line.rotate_about_origin, -t,
            ApplyFunction(rot, labela),
            ApplyFunction(rot, labelb),
            ApplyFunction(rot, labelc),
        )
        self.wait()
        text_scale = 0.625
        yz = Line(ORIGIN, ro * RIGHT).set_color(ORANGE)
        wx = Line(ro * LEFT, cid.get_center()).set_color(ORANGE)
        wd = Dot().scale(0.625).move_to(ro * LEFT).set_color(YELLOW)
        zd = Dot().scale(0.625).move_to(ro * RIGHT).set_color(YELLOW)
        xd = Dot().scale(0.625).move_to(ci.get_center() + ri * LEFT).set_color(YELLOW)
        yd = Dot().scale(0.625).move_to(ci.get_center() + ri * RIGHT).set_color(YELLOW)
        #labelo = TexMobject("O").scale(text_scale).next_to(cod.get_center(), direction = DL, buff = 0.1)
        labelw = TexMobject("W").scale(text_scale).next_to(wd.get_center(), direction = DL, buff = 0.1)
        labely = TexMobject("Y").scale(text_scale).next_to(yd.get_center(), direction = DOWN, buff = 0.1)
        labelx = TexMobject("X").scale(text_scale).next_to(xd.get_center(), direction = DOWN, buff = 0.1)
        labelz = TexMobject("Z").scale(text_scale).next_to(zd.get_center(), direction = DR, buff = 0.1)
        self.play(
            *[Write(obj) for obj in [yz, wx]],
            *[FadeInFromLarge(obj, 5) for obj in [wd, xd, yd, zd]],
        )
        self.wait()
        self.play(
            *[Write(i) for i in [labelw, labely, labelx, labelz]],
            *[FadeOut(obj) for obj in [labela, labelb, labelc]]
        )
        self.wait()
        self.play(
            self.camera_frame.shift, 2.5 * RIGHT,
            #titl.shift, 2.5 * RIGHT,
        )
        derivcr = VGroup(
            TexMobject("=\\frac{WX}{WY}\\cdot\\frac{YZ}{XZ}").scale(text_scale),
            TexMobject("=\\frac{R-d-r}{R-d+r}\\cdot\\frac{R+d-r}{R+d+r}").scale(text_scale),
            TexMobject("=\\frac{(R-r)^2-d^2}{(R+r)^2-d^2}").scale(text_scale),
        )
        derivcr[1][0][1].set_color(self.excircle_color)
        derivcr[1][0][3].set_color(BLUE)
        derivcr[1][0][5].set_color(self.incircle_color)
        derivcr[1][0][7].set_color(self.excircle_color)
        derivcr[1][0][9].set_color(BLUE)
        derivcr[1][0][11].set_color(self.incircle_color)
        derivcr[1][0][13].set_color(self.excircle_color)
        derivcr[1][0][15].set_color(BLUE)
        derivcr[1][0][17].set_color(self.incircle_color)
        derivcr[1][0][19].set_color(self.excircle_color)
        derivcr[1][0][21].set_color(BLUE)
        derivcr[1][0][23].set_color(self.incircle_color)
        derivcr[2][0][2].set_color(self.excircle_color)
        derivcr[2][0][4].set_color(self.incircle_color)
        derivcr[2][0][8].set_color(BLUE)
        derivcr[2][0][12].set_color(self.excircle_color)
        derivcr[2][0][14].set_color(self.incircle_color)
        derivcr[2][0][18].set_color(BLUE)
        lam = TexMobject("\\lambda").scale(text_scale)
        defns = VGroup(
            TexMobject("\\vert WZ \\vert =2R, \\text{ } \\vert XY \\vert =2r\\text{ and}").scale(text_scale),
            TexMobject("d\\text{ is the distance between the centers}").scale(text_scale),
        ).arrange(direction = DOWN, buff = 0.2)
        defns[0][0][6].set_color(self.excircle_color)
        defns[0][0][14].set_color(self.incircle_color)
        defns[1][0][0].set_color(BLUE)
        defns[1].align_to(defns[0], direction = LEFT)
        defns.move_to(6.5 * RIGHT + 0.75 * ro * UP)
        derivcr.arrange(direction = DOWN, buff = 0.2)
        for text in derivcr[1:]:
            text.align_to(derivcr[0], direction = LEFT)
        derivcr.move_to(1.90625 * ro * RIGHT + 0 * ro * DOWN)
        lam.next_to(derivcr[0], direction = LEFT, buff = 0.1)
        self.play(Write(defns))
        self.wait()
        self.play(Write(lam))
        for text in derivcr:
            self.play(Write(text))
            self.wait()
        self.play(
            FadeOut(defns),
            FadeOut(derivcr[0]),
            FadeOut(derivcr[1]),
            derivcr[2].next_to, lam, {"direction": RIGHT, "buff": 0.1}
        )
        fcr = VGroup(lam, derivcr[2])
        rect = SurroundingRectangle(fcr)
        self.play(ShowCreation(rect))
        fcr.add(rect)
        self.play(fcr.shift, 1.5 * RIGHT + 2.5 * UP)
        self.wait()
        fade_gp = VGroup(ca, cb, cc, labelw, labelx, labely, labelz, wd, xd, yd, zd)
        self.play(
            self.camera_frame.shift, 2 * LEFT,
            #titl.shift, 2 * LEFT,
            fade_gp.fade, fade_val,
        )
        self.wait()
        ortho = Dot().scale(0.5).move_to(RadicalAxis(ci, co).get_start())
        ext_line = DashedLine(2.5 * ortho.get_center(), wd.get_center()).set_color(ORANGE).fade(fade_val)
        tci = CircleTangentLines(ci, ortho).set_color(GREY)
        tco = CircleTangentLines(co, ortho).set_color(GREY)
        tlen = tci[0].get_length()
        orcir = Circle(radius = tlen, arc_center = ortho.get_center()).set_color(DARK_BROWN)
        #orcir = Arc(radius = tlen, arc_center = ortho.get_center(), start_angle = PI / 5, angle = -2 * PI / 5).set_color(DARK_BROWN)
        self.play(FadeInFromLarge(ortho, 5), Write(ext_line))
        self.wait()
        self.play(ShowCreation(tci), ShowCreation(tco))
        self.wait()
        self.play(
            self.camera_frame.shift, ortho.get_center(),
            self.camera_frame.scale, 1.25,
            tci.fade, 0.5,
            Write(orcir),
            tco.fade, 0.5,
        )
        self.wait()
        def perp_gp(ob):
            obj = ob.copy()
            oblen = obj.get_length()
            res_gp = VGroup()
            t_line = DashedLine(obj.get_end(), obj.get_end() + 0.75 * (obj.get_end() - obj.get_start()) / oblen).set_color(GREY)
            tt_line = t_line.copy().rotate(PI / 2, about_point = obj.get_end())
            res_gp.add(t_line)
            res_gp.add(tt_line)
            t_el = Angle(tt_line.get_end(), obj.get_end(), t_line.get_end(), radius = 0.25)
            res_gp.add(t_el)
            return res_gp
        po_gp = perp_gp(tco[1])
        pi_gp = perp_gp(tci[1])
        self.play(Write(po_gp), Write(pi_gp))
        self.wait()
        invcen = Dot().scale(0.5).move_to(ortho.get_center() + tlen * LEFT).set_color(YELLOW)
        #inve_rad = Line(invcen.get_center(), co.get_center() + ro * LEFT).get_length()
        inve_rad = ro
        invcirc = DashedVMobject(Circle(radius = inve_rad, arc_center = invcen.get_center()), num_dashes = 150).set_color(GREEN)
        self.play(FadeInFromLarge(invcen, 5), FadeOut(tci), FadeOut(tco), FadeOut(po_gp), FadeOut(pi_gp), FadeOut(ortho))
        self.wait()
        self.play(Write(invcirc))
        self.wait()
        ort_gp = perp_gp(Line(LEFT +  tlen * LEFT + ortho.get_center(), tlen * LEFT + ortho.get_center()))
        self.play(Write(ort_gp))
        self.wait()
        iort = Circle_Inversion(orcir, inversion_radius = inve_rad, about_point = invcirc.get_center()).set_color(DARK_BROWN)
        self.play(Write(iort), FadeOut(ort_gp), orcir.fade, fade_val)
        self.wait()
        ico = Circle_Inversion(co, inversion_radius = inve_rad, about_point = invcirc.get_center()).set_color(YELLOW)
        coort = perp_gp(Line(LEFT + ro * LEFT, ro * LEFT))
        self.play(FadeIn(coort), FadeIn(po_gp))
        self.wait()
        self.play(Write(ico), FadeOut(coort), FadeOut(po_gp))
        self.wait()
        ici = Circle_Inversion(ci, inversion_radius = inve_rad, about_point = invcirc.get_center()).set_color(PURPLE)
        ica = Circle_Inversion(ca, inversion_radius = inve_rad, about_point = invcirc.get_center()).set_color(RED).fade(fade_val)
        icb = Circle_Inversion(cb, inversion_radius = inve_rad, about_point = invcirc.get_center()).set_color(RED).fade(fade_val)
        icc = Circle_Inversion(cc, inversion_radius = inve_rad, about_point = invcirc.get_center()).set_color(RED).fade(fade_val)
        self.play(Write(ici))
        self.wait()
        self.play(*[Write(obj) for obj in [ica, icb, icc]])
        self.wait()
        rico = np.linalg.norm(ico.get_center() - ico.get_start())
        rici = np.linalg.norm(ici.get_center() - ici.get_start())
        zval = rico / ro
        self.play(
            self.camera_frame.move_to, ico.get_center() + rico * RIGHT / 2,
            self.camera_frame.scale, zval
        )
        wdd = Dot().scale(0.625 * zval).move_to(ico.get_center() + rico * RIGHT).set_color(YELLOW)
        zdd = Dot().scale(0.625 * zval).move_to(ico.get_center() + rico * LEFT).set_color(YELLOW)
        xdd = Dot().scale(0.625 * zval).move_to(ici.get_center() + rici * RIGHT).set_color(YELLOW)
        ydd = Dot().scale(0.625 * zval).move_to(ici.get_center() + rici * LEFT).set_color(YELLOW)
        ddgp = VGroup(wdd, zdd, xdd, ydd)
        #labelo = TexMobject("O").scale(text_scale).next_to(cod.get_center(), direction = DL, buff = 0.1)
        labelwd = TexMobject("W'").scale(text_scale * zval).next_to(wdd.get_center(), direction = DR, buff = 0.1 * zval)
        labelyd = TexMobject("Y'").scale(text_scale * zval).next_to(ydd.get_center(), direction = DOWN, buff = 0.1 * zval)
        labelxd = TexMobject("X'").scale(text_scale * zval).next_to(xdd.get_center(), direction = DOWN, buff = 0.1 * zval)
        labelzd = TexMobject("Z'").scale(text_scale * zval).next_to(zdd.get_center(), direction = DL, buff = 0.1 * zval)
        labeldd = VGroup(labelwd, labelyd, labelxd, labelzd)
        self.play(FadeInFromLarge(ddgp, 5))
        self.wait()
        self.play(Write(labeldd))
        self.wait()
        ztext_scale = text_scale * zval * 1.125
        idefns = TexMobject("\\text{Let }\\vert W'Z' \\vert =2R' \\text{ and } \\vert X'Y' \\vert =2r'").scale(ztext_scale).move_to(ico.get_center() + rico * UP + 2 * rico * RIGHT)
        idefns[0][11:13].set_color(self.excircle_color)
        idefns[0][24:].set_color(self.incircle_color)
        crderiv = VGroup(
            TexMobject("=\\frac{W'X'}{W'Y'}\\cdot\\frac{Y'Z'}{X'Z'}").scale(ztext_scale),
            TexMobject("=\\frac{R'-r'}{R'+r'}\\cdot\\frac{R'-r'}{R'+r'}").scale(ztext_scale),
            TexMobject("=\\left(\\frac{R'-r'}{R'+r'}\\right)^2").scale(ztext_scale),
            TexMobject("=\\frac{3}{4}").scale(ztext_scale)
        ).arrange(direction = DOWN, buff = 0.2 * zval)
        crderiv[1][0][1:3].set_color(self.excircle_color)
        crderiv[1][0][4:6].set_color(self.incircle_color)
        crderiv[1][0][7:9].set_color(self.excircle_color)
        crderiv[1][0][10:12].set_color(self.incircle_color)
        crderiv[1][0][13:15].set_color(self.excircle_color)
        crderiv[1][0][16:18].set_color(self.incircle_color)
        crderiv[1][0][19:21].set_color(self.excircle_color)
        crderiv[1][0][22:24].set_color(self.incircle_color)
        crderiv[2][0][2:4].set_color(self.excircle_color)
        crderiv[2][0][5:7].set_color(self.incircle_color)
        crderiv[2][0][8:10].set_color(self.excircle_color)
        crderiv[2][0][11:13].set_color(self.incircle_color)
        for text in crderiv[1:]:
            text.align_to(crderiv[0], direction = LEFT)
        crderiv.move_to(ico.get_center() + 1.9375 * rico * RIGHT + 0.25 * rico * UP)
        ilam = TexMobject("\\lambda").scale(ztext_scale).next_to(crderiv[0], direction = LEFT, buff = 0.1 * zval)
        self.play(Write(idefns))
        self.wait()
        self.play(Write(ilam))
        for text in crderiv:
            self.play(Write(text))
            self.wait()
        self.play(
            *[FadeOut(obj) for obj in [invcirc, orcir, ext_line, ortho]],
            FadeOut(crderiv[:-1]),
            crderiv[-1].next_to, ilam, {"direction": RIGHT, "buff": 0.1 * zval}
        )
        self.wait()
        self.play(
            self.camera_frame.move_to, 2 * RIGHT,
            self.camera_frame.scale, 1 / zval / 1.25,
            *[FadeOut(obj) for obj in [idefns, ilam, crderiv[-1]]]
        )
        self.wait()
        self.play(
            fcr.move_to, 6 * RIGHT + ro * UP / 2,
        )
        lamtemp = TexMobject("\\frac{3}{4}").scale(text_scale).move_to(fcr[0].get_center())
        self.play(FadeOut(rect), Transform(fcr[0], lamtemp))
        self.wait()
        res = TexMobject("d^2=R^2+r^2-14Rr", tex_to_color_map = {"r": self.incircle_color, "R": self.excircle_color, "d": BLUE}).scale(text_scale).move_to(fcr.get_center())
        '''res[0][0].set_color(BLUE)
        res[0][3].set_color(self.excircle_color)
        res[0][6].set_color(self.incircle_color)
        res[0][11].set_color(self.excircle_color)
        res[0][12].set_color(self.incircle_color)'''
        res.shift(DOWN)
        self.play(Write(res))
        self.wait(2)

class CrossRatioScenev(GraphScene, MovingCameraScene):
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
        titl = TexMobject("\\text{\\underline{An important Invariant..}}").move_to(3.5 * UP)
        inv_rad = 2
        c = Circle(radius = inv_rad).set_color(GREEN)
        od = Dot().scale(0.5)
        labelo = TexMobject("O").scale(2 / 3).next_to(od, direction = DOWN, buff = 0.1)
        ad = Dot().scale(0.5).move_to(1.25 * RIGHT)
        labela = TexMobject("A").scale(2 / 3).next_to(ad, direction = DOWN, buff = 0.1)
        bd = Dot().scale(0.5).move_to(1.5 * (np.cos(60 * DEGREES) * RIGHT + np.sin(60 * DEGREES) * UP) / 1)
        labelb = TexMobject("B").scale(2 / 3).next_to(bd, direction = RIGHT, buff = 0.1)
        dad = Circle_Inversion(ad, inversion_radius = inv_rad).scale(2 / 3)
        dbd = Circle_Inversion(bd, inversion_radius = inv_rad).scale(2 / 3)
        labelad = TexMobject("A'").scale(2 / 3).next_to(dad, direction = DOWN, buff = 0.1)
        labelbd = TexMobject("B'").scale(2 / 3).next_to(dbd, direction = RIGHT, buff = 0.1)
        loa = Line(ORIGIN, ad.get_center())
        lob = Line(ORIGIN, bd.get_center())
        load = Line(ORIGIN, dad.get_center())
        lobd = Line(ORIGIN, dbd.get_center())
        lab = DashedLine(ad.get_center(), bd.get_center()).set_color(GREY)
        ladbd = DashedLine(dad.get_center(), dbd.get_center()).set_color(GRAY)
        angoab = Angle(ORIGIN, ad.get_center(), bd.get_center(), color = RED)
        angoba = Angle(ORIGIN, bd.get_center(), ad.get_center(), color = BLUE)
        angoadbd = Angle(ORIGIN, dad.get_center(), dbd.get_center(), color = BLUE)
        angobdad = Angle(ORIGIN, dbd.get_center(), dad.get_center(), color = RED)
        self.play(Write(c), Write(od), Write(labelo), Write(titl))
        self.wait()
        ivgp = VGroup(ad, bd, dad, dbd, labela, labelb, labelad, labelbd, loa, lob, load, lobd, lab, ladbd, angoab, angoba, angoadbd, angobdad)
        self.play(FadeIn(ivgp))
        def tri_updater(obj):
            tad, tbd, tdad, tdbd, tlabela, tlabelb, tlabelad, tlabelbd, tloa, tlob, tload, tlobd, tlab, tladbd, tangoab, tangoba, tangoadbd, tangobdad = obj
            tdad.become(Circle_Inversion(tad, inversion_radius = inv_rad).scale(2 / 3))
            tdbd.become(Circle_Inversion(tbd, inversion_radius = inv_rad).scale(2 / 3))
            tlabela.next_to(tad, direction = DOWN, buff = 0.1)
            tlabelb.next_to(tbd, direction = RIGHT, buff = 0.1)
            tlabelad.next_to(tdad, direction = DOWN, buff = 0.1)
            tlabelbd.next_to(tdbd, direction = RIGHT, buff = 0.1)
            tloa.become(Line(ORIGIN, tad.get_center()))
            tlob.become(Line(ORIGIN, tbd.get_center()))
            tload.become(Line(ORIGIN, tdad.get_center()))
            tlobd.become(Line(ORIGIN, tdbd.get_center()))
            tlab.become(DashedLine(tad.get_center(), tbd.get_center()).set_color(GREY))
            tladbd.become(DashedLine(tdad.get_center(), tdbd.get_center()).set_color(GREY))
            tangoab.become(Angle(ORIGIN, tad.get_center(), tbd.get_center(), color = RED))
            tangoba.become(Angle(ORIGIN, tbd.get_center(), tad.get_center(), color = BLUE))
            tangoadbd.become(Angle(ORIGIN, tdad.get_center(), tdbd.get_center(), color = BLUE))
            tangobdad.become(Angle(ORIGIN, tdbd.get_center(), tdad.get_center(), color = RED))
        '''ivgp.add_updater(tri_updater)
        self.add(ivgp)
        self.play(ad.shift, RIGHT, run_time = 5, rate_func = there_and_back)
        ivgp.clear_updaters()'''
        self.play(
            self.camera_frame.shift, 4 * RIGHT + 0.5 * UP,
            titl.shift, 4 * RIGHT + 0.5 * UP,
        )
        text_scale = 0.625
        deriv = VGroup(
            TexMobject("\\text{By definition, }OA \\cdot OA' = OB \\cdot OB' = r^2").scale(text_scale),
            TexMobject("\\therefore \\triangle{OA'B'} \\equiv \\triangle{OBA}").scale(text_scale),
            TexMobject("\\frac{A'B'}{OB'}=\\frac{BA}{OA}\\text{ (by similarity)}").scale(text_scale),
            TexMobject("A'B'=BA \\cdot \\frac{r^2}{OA \\cdot OB}").scale(text_scale),
        ).arrange(direction = DOWN, buff = 0.5)
        deriv.shift(7 * RIGHT + 0.5 * UP)
        for text in deriv[1:]:
            text.align_to(deriv[0], direction = LEFT)
        for text in deriv:
            self.play(Write(text))
            self.wait()
        rect = SurroundingRectangle(deriv[-1])
        self.play(ShowCreation(rect))
        self.wait()
        self.play(
            self.camera_frame.shift, LEFT,
            titl.shift, LEFT,
            FadeOut(deriv[:3]),
            deriv[-1].move_to, 8 * RIGHT + 4 * UP,
            rect.move_to, 8 * RIGHT + 4 * UP,
            FadeOut(ivgp)
        )
        wd = Dot().scale(0.5).move_to(0.5 * RIGHT + 1 * UP)
        xd = Dot().scale(0.5).move_to(2.5 * LEFT + 1.625 * UP)
        yd = Dot().scale(0.5).move_to(1.5 * LEFT + 3 * UP)
        zd = Dot().scale(0.5).move_to(1.5 * RIGHT + 2.25 * UP)
        labelw = TexMobject("W").scale(text_scale).next_to(wd, direction = RIGHT, buff = 0.1)
        labely = TexMobject("Y").scale(text_scale).next_to(yd, direction = RIGHT, buff = 0.1)
        labelx = TexMobject("X").scale(text_scale).next_to(xd, direction = RIGHT, buff = 0.1)
        labelz = TexMobject("Z").scale(text_scale).next_to(zd, direction = RIGHT, buff = 0.1)
        d_gp = VGroup(wd, xd, yd, zd, labelw, labely, labelx, labelz)
        dwd = Circle_Inversion(wd, inversion_radius = inv_rad).scale(0.5)
        dxd = Circle_Inversion(xd, inversion_radius = inv_rad).scale(0.5)
        dyd = Circle_Inversion(yd, inversion_radius = inv_rad).scale(0.5)
        dzd = Circle_Inversion(zd, inversion_radius = inv_rad).scale(0.5)
        labelwd = TexMobject("W'").scale(text_scale).next_to(dwd, direction = RIGHT, buff = 0.1)
        labelyd = TexMobject("Y'").scale(text_scale).next_to(dyd, direction = RIGHT, buff = 0.1)
        labelxd = TexMobject("X'").scale(text_scale).next_to(dxd, direction = RIGHT, buff = 0.1)
        labelzd = TexMobject("Z'").scale(text_scale).next_to(dzd, direction = RIGHT, buff = 0.1)
        inv_gp = VGroup(dwd, dxd, dyd, dzd, labelwd, labelyd, labelxd, labelzd)
        self.play(Write(d_gp))
        self.wait()
        text_scale = 0.625
        cr_ra = TexMobject("\\text{Cross Ratio}").scale(text_scale).move_to(4 * RIGHT + 0.5 * UP)
        eq_sn = TexMobject("=").scale(text_scale).next_to(cr_ra, direction = RIGHT, buff = 0.1)
        defn = TexMobject("\\frac{WX}{WY}\\cdot\\frac{YZ}{XZ}").scale(text_scale).next_to(eq_sn, direction = RIGHT, buff = 0.1)
        lam = TexMobject("\\lambda").scale(text_scale).next_to(eq_sn, direction = LEFT, buff = 0.1)
        equ_grp = VGroup(cr_ra, eq_sn, defn)
        self.play(Write(equ_grp))
        self.wait()
        self.play(Transform(cr_ra, lam))
        self.wait()
        fade_val = 0.75
        self.play(
            d_gp.fade, fade_val,
            Write(inv_gp),
            equ_grp.shift, 1.5 * LEFT,
        )
        defnd = TexMobject("=\\frac{W'X'}{W'Y'}\\cdot\\frac{Y'Z'}{X'Z'}").scale(text_scale).next_to(defn, direction = RIGHT, buff = 0.1)
        self.play(Write(defnd))
        self.wait()
        equ_grp.add(defnd)
        self.play(equ_grp.shift, 1.5 * UP)
        proofpara = VGroup(
            TexMobject("\\text{\\underline{Proof:}}").scale(text_scale),
            TexMobject("\\frac{W'X'}{W'Y'}\\cdot\\frac{Y'Z'}{X'Z'}=\\frac{WX\\cdot\\frac{r^2}{OW\\cdot OX}}{WY\\cdot\\frac{r^2}{OW\\cdot OY}}\\cdot\\frac{YZ\\cdot\\frac{r^2}{OY\\cdot OZ}}{XZ\\cdot\\frac{r^2}{OX\\cdot OZ}}").scale(text_scale),
            TexMobject("=\\frac{WX}{WY}\\cdot\\frac{YZ}{XZ}").scale(text_scale),
        )
        proofpara[0].next_to(equ_grp, direction = DOWN, buff = 0.875)
        for text in proofpara:
            text.align_to(cr_ra, direction = LEFT)
        proofpara[1].shift(0.2 * DOWN)
        proofpara[2].shift(1.3 * DOWN + 1.9 * RIGHT)
        self.play(Write(proofpara))
        self.wait()
        self.play(
            self.camera_frame.shift, 3 * LEFT + 0.5 * DOWN,
            titl.shift, 3 * LEFT + 0.5 * DOWN,
            *[FadeOut(obj) for obj in [proofpara, d_gp, inv_gp, od, labelo, c, deriv[-1], rect]],
            equ_grp.move_to, ORIGIN,
            #equ_grp.shift, -e,
            equ_grp.scale, 1.5,
        )
        self.wait(2)

class InversionProperties(GraphScene, MovingCameraScene):
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
        scale_val = 0.75
        #fade_val = 0.96875
        fade_val = 0.875
        reverse_fade = 1 - 1 / (1 - fade_val)
        ttl = TextMobject("$\\text{\\underline{Properties of Inversion}}$").move_to(3.5 * UP)
        #brect = BackgroundColoredVMobjectDisplayer
        brect = BackgroundRectangle(ttl)
        #ttl.add_to_back(brect)
        titl = VGroup(brect, ttl)
        self.add_foreground_mobject(titl)
        prop_list = VGroup(
            TextMobject("Inverse of a point $A$ is defined such that $OA\\cdot OA'=r^2$", tex_to_color_map = {"$OA\\cdot OA'=r^2$": YELLOW}).scale(scale_val),
            TextMobject("Lines passing through $O$ remains unchanged under Inversion").scale(scale_val),
            TextMobject("Circles passing through $O$ inverts to Lines not passing through $O$").scale(scale_val),
            TextMobject("Circles not passing through $O$ inverts to Circles not passing through $O$").scale(scale_val),
            TextMobject("Circles orthogonal to Inversion circle remain unchanged").scale(scale_val),
            TextMobject("Angle of Intersection between two curves remain unchanged under Inversion").scale(scale_val),
        ).arrange(direction = DOWN, buff = 0.6875)
        for text in prop_list[1:]:
            text.align_to(prop_list[0], direction = LEFT)
        prop_list.to_edge()
        prop_list.shift(0.5 * DOWN)
        self.play(Write(titl))
        self.wait()
        inv_rad = 2
        c = Circle(radius = inv_rad).set_color(GREEN)
        od = Dot().scale(0.5)
        labelo = TexMobject("O").scale(2 / 3).next_to(od, direction = DOWN, buff = 0.1)
        ad = Dot().scale(0.5).move_to(1.25 * RIGHT)
        labela = TexMobject("A").scale(2 / 3).next_to(ad, direction = DOWN, buff = 0.1)
        bd = Dot().scale(0.5).move_to(1.5 * (np.cos(60 * DEGREES) * RIGHT + np.sin(60 * DEGREES) * UP) / 1)
        labelb = TexMobject("B").scale(2 / 3).next_to(bd, direction = RIGHT, buff = 0.1)
        dad = Circle_Inversion(ad, inversion_radius = inv_rad).scale(2 / 3)
        dbd = Circle_Inversion(bd, inversion_radius = inv_rad).scale(2 / 3)
        labelad = TexMobject("A'").scale(2 / 3).next_to(dad, direction = DOWN, buff = 0.1)
        labelbd = TexMobject("B'").scale(2 / 3).next_to(dbd, direction = RIGHT, buff = 0.1)
        rayoa = Line(ORIGIN, 10 * ad.get_center()).set_color(GREY)
        loa = Line(ORIGIN, ad.get_center())
        lob = Line(ORIGIN, bd.get_center())
        load = Line(ORIGIN, dad.get_center())
        lobd = Line(ORIGIN, dbd.get_center())
        lab = DashedLine(ad.get_center(), bd.get_center()).set_color(GREY)
        ladbd = DashedLine(dad.get_center(), dbd.get_center()).set_color(GRAY)
        angoab = Angle(ORIGIN, ad.get_center(), bd.get_center(), color = RED)
        angoba = Angle(ORIGIN, bd.get_center(), ad.get_center(), color = BLUE)
        angoadbd = Angle(ORIGIN, dad.get_center(), dbd.get_center(), color = BLUE)
        angobdad = Angle(ORIGIN, dbd.get_center(), dad.get_center(), color = RED)
        self.play(Write(c), Write(od), Write(labelo))
        self.wait()
        origgp = VGroup(c, od, labelo)
        self.play(
            origgp.fade, fade_val,
            Write(prop_list[0])
        )
        self.wait()
        self.play(
            prop_list[0].fade, fade_val,
            FadeIn(ad), FadeIn(labela), FadeIn(rayoa),
            origgp.fade, reverse_fade,
        )
        self.wait()
        self.play(
            FadeIn(dad), FadeIn(labelad),
        )
        apath = VMobject()
        apath.set_points_smoothly([ad.get_center(), 2 * RIGHT + 1 * UP, 3 * RIGHT, ad.get_center()])
        #apath.set_points_smoothly([ad.get_center(), 2 * RIGHT + 1 * UP, 2.1 * UP, 0.5 * LEFT + 2 * DOWN, 2.5 * LEFT, ad.get_center()])
        apath.fade(0.99)
        #bpath = RoseCurve(radius = 2, kval = 3 / 2)
        self.add(apath)
        a_gp = VGroup(ad, labela, dad, labelad, rayoa)
        def pt_upd(obj):
            if len(obj) == 5:
                pt, ptl, pti, ptil, raypt = obj
            else:
                pt, ptl, pti, ptil = obj
            ptl.next_to(pt, direction = DOWN, buff = 0.1)
            pti.become(Circle_Inversion(pt, inversion_radius = inv_rad).scale(2 / 3))
            ptil.next_to(pti, direction = DOWN, buff = 0.1)
            if len(obj) == 5:
                raypt.become(Line(ORIGIN, 10 * pt.get_center()).set_color(GREY))
        a_gp.add_updater(pt_upd)
        self.add(a_gp)
        self.play(MoveAlongPath(ad, apath), run_time = 5, rate_func = linear)
        a_gp.clear_updaters()
        self.wait()
        self.play(
            *[FadeOut(obj) for obj in [ad, dad, labela, labelad, rayoa]],
            origgp.fade, fade_val,
        )
        self.wait()
        a_gp.remove(rayoa)
        self.play(FadeIn(prop_list[1]),)
        self.wait()
        self.play(
            prop_list[1].fade, fade_val,
            origgp.fade, reverse_fade,
        )
        self.wait()
        exline = Line(20 * RIGHT + 10 * UP, 20 * LEFT + 10 * DOWN)
        self.play(Write(exline))
        self.wait()
        ad.move_to(4 * LEFT + 2 * DOWN)
        dad.become(Circle_Inversion(ad, inversion_radius = inv_rad).scale(2 / 3))
        labela.next_to(ad, direction = DOWN, buff = 0.1)
        labelad.next_to(dad, direction = DOWN, buff = 0.1)
        a_gp.add_updater(pt_upd)
        self.play(
            *[FadeIn(obj) for obj in [ad, dad, labela, labelad]],
        )
        self.add(a_gp)
        self.play(ad.move_to, 1 * LEFT + 0.5 * DOWN, run_time = 5, rate_func = smooth)
        a_gp.clear_updaters()
        self.wait()
        self.play(
            FadeOut(a_gp),
            FadeOut(exline),
            origgp.fade, fade_val,
        )
        self.wait()
        self.play(FadeIn(prop_list[2]))
        self.wait()
        self.play(
            origgp.fade, reverse_fade,
            prop_list[2].fade, fade_val,
        )
        self.wait()
        ccir = Circle(radius = 0.8).shift(0.8 * RIGHT)
        cciri = Circle_Inversion(ccir, inversion_radius = inv_rad)
        self.play(Write(ccir), Write(cciri))
        self.wait()
        ad.move_to(ccir.get_start())
        dad.become(Circle_Inversion(ad, inversion_radius = inv_rad).scale(2 / 3))
        labela.next_to(ad, direction = DOWN, buff = 0.1)
        labelad.next_to(dad, direction = DOWN, buff = 0.1)
        a_gp.add_updater(pt_upd)
        self.play(
            FadeIn(a_gp),
        )
        self.add(a_gp)
        self.play(MoveAlongPath(ad, Arc(radius = 0.8, arc_center = ccir.get_center())), run_time = 5, rate_func = there_and_back)
        a_gp.clear_updaters()
        self.wait()
        self.play(
            a_gp.fade, fade_val,
            ccir.fade, fade_val,
            cciri.fade, fade_val,
            origgp.fade, fade_val,
        )
        self.wait()
        self.play(FadeIn(prop_list[3]))
        self.wait()
        self.play(
            origgp.fade, reverse_fade,
            a_gp.fade, reverse_fade,
            ccir.fade, reverse_fade,
            cciri.fade, reverse_fade,
            prop_list[3].fade, fade_val,
        )
        self.wait()
        in_gp = VGroup(ccir, cciri, ad, labela, dad, labelad)
        def gp_upd(obj):
            tcir, tciri, tad, tla, tdad, tlad = obj
            tciri.become(Circle_Inversion(tcir, inversion_radius = inv_rad))
            tad.become(Dot().scale(0.5).move_to(tcir.get_start()))
            tdad.become(Circle_Inversion(tad, inversion_radius = inv_rad).scale(2 / 3))
            tla.next_to(tad, direction = DOWN, buff = 0.1)
            tlad.next_to(tdad, direction = DOWN, buff = 0.1)
        in_gp.add_updater(gp_upd)
        self.add(in_gp)
        self.play(
            ccir.shift, 1.6 * RIGHT,
            ccir.rotate_about_origin, PI / 3,
            run_time = 5
        )
        self.wait()
        in_gp.clear_updaters()
        self.play(
            in_gp.fade, fade_val,
            origgp.fade, fade_val,
        )
        self.wait()
        self.play(FadeIn(prop_list[4]))
        self.wait()
        self.play(
            in_gp.fade, reverse_fade,
            origgp.fade, reverse_fade,
            prop_list[4].fade, fade_val,
        )
        self.wait()
        ortpt = Dot().scale(0.5).move_to(4 * RIGHT + 2 * DOWN)
        tanlin = CircleTangentLines(c, ortpt)
        ortan = VGroup()
        for obj in tanlin:
            ortan.add(DashedVMobject(obj, num_dashes = 45))
        ortan.set_color(GREY)
        ccen = Dot().scale(0.5).move_to(ccir.get_center())
        self.play(Write(ortpt), Write(ortan), Write(ccen))
        self.wait()
        ccir.move_to(ccen.get_center())
        in_gp.add(ccen)
        def gpc_upd(obj):
            tcir, tciri, tad, tla, tdad, tlad, tcen = obj
            tcir.become(Circle(radius = np.linalg.norm(tcen.get_center() - tad.get_center())).move_to(tcen.get_center()))
            tciri.become(Circle_Inversion(tcir, inversion_radius = inv_rad))
            #tad.become(Dot().scale(0.5).move_to(tcir.get_start()))
            tdad.become(Circle_Inversion(tad, inversion_radius = inv_rad).scale(2 / 3))
            tla.next_to(tad, direction = DOWN, buff = 0.1)
            tlad.next_to(tdad, direction = DOWN, buff = 0.1)
        in_gp.add_updater(gpc_upd)
        self.add(in_gp)
        self.play(ccen.move_to, ortpt.get_center(), ad.move_to, tanlin[0].get_end(), run_time = 5)
        in_gp.clear_updaters()
        self.wait()
        def perp_gp(ob):
            obj = ob.copy()
            oblen = obj.get_length()
            res_gp = VGroup()
            t_line = DashedLine(obj.get_end(), obj.get_end() + 0.75 * (obj.get_end() - obj.get_start()) / oblen).set_color(GREY)
            tt_line = t_line.copy().rotate(PI / 2, about_point = obj.get_end())
            res_gp.add(t_line)
            res_gp.add(tt_line)
            t_el = Angle(tt_line.get_end(), obj.get_end(), t_line.get_end(), radius = 0.25)
            res_gp.add(t_el)
            return res_gp
        pgp = perp_gp(Line(ORIGIN, ad.get_center()))
        self.play(Write(pgp))
        self.wait()
        self.play(
            in_gp.fade, fade_val,
            origgp.fade, fade_val,
            *[FadeOut(obj) for obj in [ortan, ccen, pgp, ortpt]],
        )
        self.wait()
        self.play(FadeIn(prop_list[5]))
        self.wait()
        othercircle = Circle().move_to(3 * UR)
        otherinv = Circle_Inversion(othercircle, inversion_radius = inv_rad)
        self.play(
            ccir.fade, reverse_fade,
            cciri.fade, reverse_fade,
            origgp.fade, reverse_fade,
            *[FadeOut(obj) for obj in [labela, labelad, ad, dad]],
            prop_list[5].fade, fade_val,
            *[Write(obj) for obj in [othercircle, otherinv]]
        )
        self.wait()
        circlesgp = VGroup(ccir, cciri, othercircle, otherinv)
        def new_upd(obj):
            c1, i1, c2, i2 = obj
            i1.become(Circle_Inversion(c1, inversion_radius = inv_rad))
            i2.become(Circle_Inversion(c2, inversion_radius = inv_rad))
        circlesgp.add_updater(new_upd)
        self.add(circlesgp)
        self.play(
            ccir.shift, 2.5 * UP + 0.75 * LEFT,
            ccir.scale_in_place, 1 / 4,
            othercircle.shift, 1.5 * DOWN + 0.25 * LEFT,
            othercircle.scale_in_place, 0.875,
            run_time = 5
        )
        circlesgp.clear_updaters()
        interpts = Intersection(ccir, othercircle)[0]
        templine = Line(ccir.get_center(), interpts.get_center())
        templine.rotate(1e-3, about_point = ccir.get_center())
        fline = Line(interpts.get_center(), interpts.get_center() - 500 * (templine.get_end() - interpts.get_center())).set_color(GREY)
        templine = Line(othercircle.get_center(), interpts.get_center())
        templine.rotate(1e-3, about_point = othercircle.get_center())
        sline = Line(interpts.get_center(), interpts.get_center() + 500 * (templine.get_end() - interpts.get_center())).set_color(GREY)
        ang = Angle(fline.get_end(), interpts.get_center(), sline.get_end(), radius = 0.25, color = YELLOW)
        self.play(Write(sline), Write(fline), Write(ang))
        self.wait()
        self.remove(cciri)
        cciri = Circle_Inversion(ccir, inversion_radius = inv_rad)
        self.add(cciri)
        interpts = Intersection(cciri, otherinv)[1]
        templine = Line(cciri.get_center(), interpts.get_center())
        templine.rotate(1e-3, about_point = cciri.get_center())
        ifline = Line(interpts.get_center(), interpts.get_center() + 1000 * (templine.get_end() - interpts.get_center())).set_color(GREY)
        templine = Line(otherinv.get_center(), interpts.get_center())
        templine.rotate(1e-3, about_point = otherinv.get_center())
        isline = Line(interpts.get_center(), interpts.get_center() - 1000 * (templine.get_end() - interpts.get_center())).set_color(GREY)
        iang = Angle(ifline.get_end(), interpts.get_center(), isline.get_end(), radius = 0.25, color = YELLOW)
        self.play(Write(isline), Write(ifline), Write(iang))
        self.wait()
        fade_list = []
        for obj in [sline, fline, isline, ifline, ang, iang, ccir, cciri, otherinv, othercircle, origgp]:
            fade_list += [obj.fade, fade_val]
        self.play(
            prop_list.fade, reverse_fade,
            *fade_list,
        )
        self.wait(2)


if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    command_B = module_name + " " + "ProblemIntroduction" + " -ps -n 20,20"
    #command_B = module_name + " " + "IntroductionScene" + " -p",
    #command_B = module_name + " -p"
    os.system(clear_cmd)
    os.system(command_A + command_B)