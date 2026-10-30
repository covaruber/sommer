(function () {
  "use strict";

  const home = document.querySelector(".sommer-home");
  if (!home) return;

  const menu = home.querySelector(".site-language");
  const buttons = menu.querySelectorAll("[data-language]");
  const translations = Array.from(home.querySelectorAll("[data-es]"), function (element) {
    return { element: element, en: element.textContent, es: element.dataset.es };
  });
  const englishTitle = document.title;
  const storageKey = "sommer-home-language";

  function applyLanguage(language) {
    translations.forEach(function (translation) {
      translation.element.textContent = translation[language];
    });
    document.documentElement.lang = language;
    document.title = language === "es"
      ? "sommer - Modelos mixtos para genética cuantitativa • sommer"
      : englishTitle;
    home.querySelector(".site-nav").setAttribute("aria-label",
      language === "es" ? "Navegación principal" : "Main navigation");
    menu.setAttribute("aria-label",
      language === "es" ? "Idioma" : "Language");
    buttons.forEach(function (button) {
      button.setAttribute("aria-pressed", String(button.dataset.language === language));
    });
  }

  let language = "en";
  try {
    if (localStorage.getItem(storageKey) === "es") language = "es";
  } catch (error) {
    console.warn("sommer: unable to read the saved homepage language.", error);
  }
  applyLanguage(language);
  menu.hidden = false;

  buttons.forEach(function (button) {
    button.addEventListener("click", function () {
      applyLanguage(button.dataset.language);
      try {
        localStorage.setItem(storageKey, button.dataset.language);
      } catch (error) {
        console.warn("sommer: unable to save the homepage language.", error);
      }
    });
  });
}());
