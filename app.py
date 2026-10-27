import streamlit as st

st.set_page_config(
    page_title="Навчальний симулятор rAAV",
    layout="wide"
)

st.title("Навчальний симулятор rAAV")

st.success("Сайт працює! Ви успішно запустили застосунок.")

st.write(
    "Це тестова сторінка. "
    "Наступним кроком додамо навчальну модель."
)

volume = st.slider(
    "Навчальний приклад: об’єм біореактора, л",
    min_value=200,
    max_value=2000,
    value=200,
    step=100
)

st.metric("Вибраний об’єм", f"{volume} л")

st.caption(
    "Цей повзунок поки лише перевіряє роботу сайту "
    "і не розраховує вихід вірусного вектора."
)
